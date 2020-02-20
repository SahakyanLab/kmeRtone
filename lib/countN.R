countN <- function(genome, base, width, ncpu=1, exclude.mito=TRUE) {
  # Calculate percentage of base N. N can be a single base e.g. {A,C,G,T} and OR set e.g. C or G (C|G) 
  # The calculation is done in a sliding window manner.
  #
  # genome    <string>      A variable name of <genome> class object.
  # base     <character>    A base e.g. {A,C,G,T}. It can be a vector as well to count a multiple occurence.
  # width     <numeric>     A width of DNA sequence to window slide.
  # ncpu      <numeric>     A number of cpu to use.
  
  # Dependencies
  #    Package   : stringi, foreach, doParallel
  #    Function  : distributeChunk
  
  #width = c(100, 500, 1000, 5000, 10000, 50000, 100000)
  
  if (exclude.mito) {
    mito = c("chrM", "chrm", "chrmt", "chrMT", "chrMt", "chrMito", "chrmito")
    chromosome.names <- names(genome)[!names(genome) %in% mito]  # exclude mitochondria
  } else {
    chromosome.names <- names(env[[genome]])
  }
  
  # total number of nucleotides inside bins to process at one time
  threshold.bases = 1e+7
  
  # initiate GC table
  N.table <- lapply(width, function(w) {
    N.table <- rep(0, 101)
    names(N.table) <- 0:100
    return(N.table)
  })
  names(N.table) <- width
  
  for (w in width) {
    
    # total bins to process at one time
    threshold.bins <- threshold.bases %/% w
    
    for (chr in chromosome.names) {
      
      len <- attr(env[[genome]], "length")[[chr]]
      
      # total bins
      total.bins <- len - w + 1
      
      chunks <- total.bins %/% threshold.bins
      
      remainder.bin <- total.bins %% threshold.bins
      
      cat("Counting N content [", base, "] of", chr, "at", w, "width...")
      start.time <- Sys.time()
      
      if (ncpu == 1) {
        
        i <- 0
        while (i < chunks) {
          
          start.bin <- i * threshold.bins + 1
          end.bin <- (i+1) * threshold.bins
          
          dna.seq <- stri_sub(genome[[chr]], start.bin:end.bin, length = w)
          
          N.count <- lapply(base, function(nt) stri_count_fixed(dna.seq, nt))
          N.count <- as.integer(Reduce(`+`, N.count) / w * 100)
          N.count <- table(N.count)
          
          N.table[[as.character(w)]][names(N.count)] <- N.table[[as.character(w)]][names(N.count)] + N.count
          
          i <- i+1
        }
        
      } else if (ncpu > 1) {
        
        cpu.chunks <- distributeChunk(chunks, ncpu)
        
        chr.start <- c(0, cpu.chunks$end * threshold.bins)+1
        chr.start <- chr.start[-length(chr.start)]
        chr.end <- cpu.chunks$end * threshold.bins
        chr.seq <- stri_sub(env[[genome]][[chr]], chr.start, chr.end+w-1)
        
        to.export <- c("base", "cpu.chunks", "threshold.bins", "w")
        
        N.tables <- foreach(n=1:length(chr.seq), chr.seq=chr.seq, .noexport=ls()[!ls() %in% to.export], .packages="stringi") %dopar% {
          
          #initiate table
          N.table <- rep(0 , 101)
          names(N.table) <- 0:100
          
          # shift index because of dna cutting
          if (n != 1) {
            cpu.chunks$end[n] <- cpu.chunks$end[n] - cpu.chunks$start[n] + 1 
          }
          
          i <- 0
          while (i < cpu.chunks$end[n]) {
            
            start.bin <- i * threshold.bins + 1
            end.bin <- (i+1) * threshold.bins
            
            dna.seq <- stri_sub(chr.seq, start.bin:end.bin, length = w)
            
            N.count <- lapply(base, function(nt) stri_count_fixed(dna.seq, nt))
            N.count <- as.integer(Reduce(`+`, N.count) / w * 100)
            
            N.count <- table(N.count)
            
            N.table[names(N.count)] <- N.table[names(N.count)] + N.count
            
            i <- i+1
          }
          return(N.table)
        }
        rm(chr.seq)
        
        N.table[[as.character(w)]] <- N.table[[as.character(w)]] + Reduce(`+`, N.tables)
      }
      
      if (remainder.bin > 0) {
        
        start.bin <- chunks * threshold.bins + 1
        end.bin <- len-w+1
        
        dna.seq <- stri_sub(env[[genome]][[chr]], start.bin:end.bin, length = w)
        
        N.count <- lapply(base, function(nt) stri_count_fixed(dna.seq, nt))
        
        N.count <- as.integer(Reduce(`+`, N.count) / w * 100)
        
        N.count <- table(N.count)
        
        N.table[[as.character(w)]][names(N.count)] <- N.table[[as.character(w)]][names(N.count)] + N.count
      }
      
      time.diff <- Sys.time() - start.time
      cat("DONE! ---", time.diff, attr(time.diff, "units"), "\n")
    }
  }
  
  if (length(N.table) == 1) N.table <- N.table[[1]]
  
  return(N.table)
}

