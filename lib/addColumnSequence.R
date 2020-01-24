addColumnSequence <- function(genomic.coordinate, genome.path, chromosome.suffix=".fa.gz", ncpu=1) {
  # Because damage sequence is expected to be 1 or 2, we can do vectorisation of 
  # substring() function and get away from memory problem.
  #
  # Because data.table is already multithreading, we can only see faster speed after 10 mil table row 
  # when using foreach loop
  
  dt <- genomic.coordinate
  
  # expand "*" to "+" and "-" if any
  if (nrow(dt[strand == "*"]) > 0) {
    dt <- rbind(dt[strand %in% c("+", "-")],
                dt[strand == "*"][, strand := "+"],
                dt[strand == "*"][, strand := "-"])
  }
  
  if (ncpu == 1) {
    
    # get sequence
    dt[, sequence := readGenome(genome.path, chromosome, start, end, form = "string"), by = chromosome]
    
    # reverse complement for minus strand
    dt[strand == "-", sequence := reverseComplement(sequence, form = "string")]
    
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
    
    dt.list <- lapply(1:length(ith), function(i) dt[ ith[i]:jth[i] ])
    
    # map the sequence
    
    dt <- foreach(dt=dt.list, .combine = "rbind", .packages = "data.table",
            .export = c("readGenome", "reverseComplement", "chromosome.suffix")) %dopar% {
      
      # get sequence
      dt[, sequence := readGenome(genome.path, chromosome, start, end, form = "string"),
         by = chromosome]
      
      # reverse complement for minus strand
      dt[strand == "-", sequence := reverseComplement(sequence, form = "string")]
      
      return(dt)
    }
    
    stopCluster(cl)
    
  }
  
  return(dt)
}