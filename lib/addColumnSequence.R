addColumnSequence <- function(genomic.coordinate, genome.path, chromosome.suffix=".fa.gz", ncpu=1) {
  
  dt <- genomic.coordinate
  
  # expand "*" to "+" and "-"
  if (nrow(dt[strand == "*"]) > 0) {
    dt <- rbind(dt[strand %in% c("+", "-")],
                dt[strand == "*"][, strand := "+"],
                dt[strand == "*"][, strand := "-"])
  }
  
  if (ncpu == 1) {
    
    dt[, sequence := lapply(1:.N, function(i) {
      
      dna.seq <- readGenome(genome.path, chromosome[i], start[i], end[i], chromosome.suffix)
      
      # if - strand, reverse complement
      if (strand[i] == "-") {
        dna.seq <- reverseComplement(dna.seq)
      }
      
      # if more than one bases, collapse to string
      if (end[i] - start[i] + 1 > 1) {
        dna.seq <- paste(dna.seq, collapse = "")
      }
      
      return(dna.seq)
      
    })]
    
  } else {
    
    require(foreach)
    require(doParallel)
    
    cl <- makeCluster(ncpu)
    registerDoParallel(cl)
    
    # distribute rows to cpu
    total.rows <- nrow(dt)
    dt.segment <- rep(total.rows %/% ncpu, ncpu)
    remainder <- total.rows %% ncpu
    if (remainder > 0) {
      remainder.rows <- rep(0, ncpu)
      remainder.rows[1:remainder] <- rep(1, remainder)
      dt.segment <- dt.segment + remainder.rows
    }
    
    ith <- c(1, cumsum(dt.segment) + 1)[1:ncpu]
    jth <- cumsum(dt.segment)
    
    dt[, sequence := NA]
    
    # map the sequence
    dt[,
       sequence := foreach(ith=ith, jth=jth, .combine = "c", .packages = "data.table",
                           .export = c("readGenome", "reverseComplement", "genome.path",
                                       "chromosome.suffix")) %dopar% {
                             
                             dna.seq <- lapply(ith:jth, function(i) {
                               
                               dna.seq <- readGenome(genome.path, chromosome[i], start[i], end[i], chromosome.suffix)
                               
                               # if - strand, reverse complement
                               if (strand[i] == "-") {
                                 dna.seq <- reverseComplement(dna.seq)
                               }
                               
                               # if more than one bases, collapse to string
                               if (end[i] - start[i] + 1 > 1) {
                                 dna.seq <- paste(dna.seq, collapse = "")
                               }
                               
                               return(dna.seq)
                             })
                             
                             return(dna.seq)
                           }
       ]
    
    stopCluster(cl)
    
  }
  
  return(dt)
}