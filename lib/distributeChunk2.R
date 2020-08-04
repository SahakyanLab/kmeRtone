distributeChunk2 <- function(numeric.vector, max.value, output) {
  # Perform a cumsum and if the max.value is reached, restart the cumsum
  #   The goal is to get distribution of group that their sum is ~ max.value
  #
  # This function was written to select a group of kmers from different region
  # that makes up to 10,000,000 kmers to be memory efficient. (only ~1GB invoked)
  # e.g. c(123, 10, 400, 200, ...) # every elements is total of kmers from a region.
  #
  # output <string>
  #    "index"   : a list of 2 element (start and end) e.g. list(start=c(1,2,3), end=c(1,2,3))
  #    "label"   : a label for the vector e.g. 1 1 1 1 1 1 1 2 2 2 2 3 3 4 4 ...
  
  # Dependency:
  #   Package: data.table - need shift() for vector fast shifting.
  
  cumsum.threshold <- cumsum(numeric.vector)
  i <- 1 # ith segment
  
  while (TRUE) {
    idx.threshold <- cumsum.threshold %/% max.value > 0
    
    if (sum(idx.threshold) == length(idx.threshold)) break() # all number bigger than max.value
    
    ith.segment <- which(idx.threshold)[i]
    if (is.na(ith.segment)) break() # NA because the the cumsum at the end cannot reach max.value
    idx.threshold[1:ith.segment] <- FALSE
    
    cumsum.threshold[idx.threshold] <- cumsum(numeric.vector[idx.threshold])
    
    # to break loop
    if (sum(idx.threshold) == 0) break() # reach at the end and cumsum reach max.value
    
    i <- i+1
  }
  
  if (output == "label") {
    return(cumsum(shift(idx.threshold, n=1, fill=FALSE))+1)
  } else if (output == "index") {
    
    if (sum(idx.threshold) == 0) {
      idx.threshold <- cumsum.threshold %/% max.value > 0
    }
    
    end <- which(idx.threshold)
    
    # if numeric(0)
    if (length(end) == 0) {
      end <- length(numeric.vector)
    }
    
    if (end[length(end)] != length(numeric.vector)) {
      end[length(end)+1] <- length(numeric.vector)
    }
    
    start <- shift(end, n=1, fill=0) + 1
    idx <- list(start = start, end = end)
    
    return(idx) 
  }
}