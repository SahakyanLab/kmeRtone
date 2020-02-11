Ncontent <- function(count.genome=F, count.damage=F, N=NULL) {
  # To calculate N percentage of damage and damage pattern in the genome
  #     at several window size. N is either G or GC
  # 1) To load pre-calculated N percentage of the genome or calculate from scratch
  # 3) To plot the density of N distribution.
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, chromosome.names, genome.path, ncpu
  #     Packages          : data.table, foreach, doParallel
  #     Functions         : readGenome, reverseComplement
  
  if (nchar(N) == 1) {
    N.bool.eval <- sprintf("dna.seq == %s", N)
  } else {
    base <- strsplit(N, "")[[1]]
    N.bool.eval <- sprintf("dna.seq == '%s' | dna.seq == '%s'", base[1], base[2])
  }
  
  width <- c(100, 500, 1000, 5000, 10000, 50000, 100000)

  
  # -------------------------------------------------------------------------------------------
  # N content of genome
  if (count.genome == TRUE) {
    
    N.table <- lapply(1:length(width), function(i) {
      N.table <- rep(0, 101)
      names(N.table) <- 0:100
      return(N.table)
    })
    names(N.table) <- width
    
    Time.Start <- Sys.time()
    for (w in width) {
      time.Start <- Sys.time()
      for (chr in chromosome.names) {
        time.start <- Sys.time()
        cat("Counting", N, "of", chr, "at", as.integer(w), "width\t")
        
        # get chromosome sequence
        dna.seq <- strsplit(genome[[chr]], "", fixed = TRUE)[[1]]
        
        # get N boolean
        N.bool <- eval(parse(text = N.bool.eval))
        
        # we don't need the sequence anymore. discard!
        rm(dna.seq)
        gc()
        
        N.count <- countSlidingBool(N.bool, w, ncpu)
        
        N.update <- c(N.table[[as.character(w)]], N.count)
        
        N.table[[as.character(w)]] <- tapply(N.update, names(N.update), sum)
        
        time.diff <- Sys.time() - time.start
        cat(time.diff[[1]], attr(time.diff, "units"), "\n")
      }
      time.Diff <- Sys.time() - time.Start
      cat("Total time taken for", w, "width is", time.Diff[[1]], attr(time.Diff, "units"), "\n\n")
    }
    Time.Diff <- Sys.time() - Time.Start
    cat("Total time taken is", Time.Diff[[1]], attr(Time.Diff, "units"), "\n")
    
    fwrite(N.table, "data/analysis/GC_genome.csv")
  }
  
  # -------------------------------------------------------------------------------------------
  # N content of the damage coordinate
  # Retrieve the pre-calculated N percentage of the genome or recalculate from scratch
  
  if (count.damage == TRUE) {
    
    genomic.coordinate[, `:=`(original_start = start, original_end = end) ]
    
    N.table <- lapply(1:length(width), function(i) {
      N.table <- rep(0, 101)
      names(N.table) <- 0:100
      return(N.table)
    })
    names(N.table) <- width
    
    for (w in width) {
      
      # shrink start and end coordinates to point to their midpoint
      genomic.coordinate[, `:=`(start = start + ((end-start+1)%/%2) - 1,
                                end = end - as.integer((end-start+1)/2 + 0.5) + 1)]
      
      # scale up to width; if middle point is odd number expand w/2 on both side; if even expand w/2-1
      genomic.coordinate[, `:=`(
        start = {
          len <- end-start+1
          s <- start
          s[len == 2] <- s - w/2 + 1
          s[len == 1] <- s - w/2
          s
        },
        end = {
          len <- end-start+1
          e <- end
          e[len == 2] <- e + w/2 - 1
          e[len == 1] <- e + w/2
          e
        }
      )]
      
      # assign NA to out of range coordinate
      trimGenCoordinate("genome.coordinate", "genome", kmertone.env)
      
      # sort by chromosome
      setkey(genomic.coordinate, chromosome)
      
      N.table <- genomic.coordinate[(!is.na(start) | (!is.na(end))), {
        
        # get chromosome sequence
        dna.seq <- strsplit(genome[[chromosome]], "")[[1]]
        
        # get N boolean
        N.bool <- eval(parse(text = N.bool.eval))
        
        # we don't need the sequence anymore. discard!
        rm(dna.seq)
        gc()
        
        N.table <- rep(0, 101)
        for (i in seq_along(start)) {
          
          N.count <- as.integer(sum(N.bool[start:end]) / (end-start+1) * 100)
          
          # update N.table
          N.table[[N.count-1]] <- N.table[[N.count-1]] + 1

        }
        
        # N.count <- mapply(function(start, end) {as.integer(sum(N.bool[start:end]) / (end-start+1) * 100)},
        #                   start, end)
        
      }, by = chromosome][[2]]
      
      # name the count
      names(N.table) <- rep(0:100, length(N.table)/101)
      
      # add same count name
      N.table[[as.character(w)]] <- tapply(N.table[[as.character(w)]], names(N.table[[as.character(w)]]), sum)
      
      # revert coordinate back to original
      genomic.coordinate[, `:=`(start = original_start, end := original_end)]
      
    }
    fwrite(N.table, "table_GC.csv")
  }
}