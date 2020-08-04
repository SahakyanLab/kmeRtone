scaleCountSliding <- function(count.vector.1, width.1, width.2, ncpu=1) {
  # This function calculates pattern matching count of sliding windows by
  # scaling up previous count of sliding window. It save times by reducing
  # the length of subsetting if we were to calculate from the beginning. 
  # Even c++ is slow for subsetting a large vector e.g. A[ 1:100000000 ]
  #      p/s: My knowledge of c++ is limited, so the last statement can be wrong.
  #           but the c++ version is slow when subsetting a large vector.
  #
  # The parallel version is not memory efficient because is has to memory copy
  # a large count vector to each cluster. 
  
  # Dependency: fcoalesce() function from data.table package.

  # check if width 2 can be scale up from width 1
  if (width.2 %% width.1 > 0) {
    stop(paste0(width.2, "can not be scaled up from", width.1))
  }
  
  total.sequence <- length(count.vector.1) + width.1 - 1
  
  if (ncpu == 1) {
    
    scale <- width.2 / width.1
    count.vector.2 <- rep(NA, total.sequence - width.2 + 1 )
    
    i <- 1
    while (i <= width.2) {
      
      # find non-overlap bin as if binning method
      bin1 <- seq(i, length(count.vector.1), width.1)
      
      # last bin for vector 2 
      bin2.idx.end <- length(bin1) %/% scale * scale

      # absorb suitable bin1 to bin2
      bin2.idx <- bin1[1:bin2.idx.end]
      
      # absorb the count
      bin2 <- count.vector.1[bin2.idx]
      
      # merge bin1
      bin2 <- colSums(matrix(bin2, nrow = scale))
      
      # position the bin appropriately to reflect sliding window
      count.vector.2[seq(i, length(count.vector.2), width.2)] <- bin2
      
      i<-i+1
    }
    
  } else {
    
    # distribute work to cpu
    chunks <- distributeChunk(width.2, ncpu)
    
    scale <- width.2 / width.1
    count.vector.2 <- rep(NA, total.sequence - width.2 + 1 )
    
    count.vector.2 <- foreach(start=chunks$start, end=chunks$end) %dopar% {
      
      i <- start
      while (i <= end) {
        
        # find non-overlap bin as if binning method
        bin1 <- seq(i, length(count.vector.1), width.1)
        
        # last bin for vector 2 
        bin2.idx.end <- length(bin1) %/% scale * scale
        
        # absorb suitable bin1 to bin2
        bin2.idx <- bin1[1:bin2.idx.end]
        
        # absorb the count
        bin2 <- count.vector.1[bin2.idx]
        
        # merge bin1
        bin2 <- colSums(matrix(bin2, nrow = scale))
        
        # position the bin appropriately to reflect sliding window
        count.vector.2[seq(i, length(count.vector.2), width.2)] <- bin2
        
        i <- i + 1
      }
            
      return(count.vector.2)
    }
    
    # merge the count
    count.vector.2 <- fcoalesce(count.vector.2)
  }
  
  
  return(count.vector.2)
}