countSlidingBool <- function(bool.vector, width, ncpu=1) {
  # This function calculate the "TRUE" percentage of sliding windows of a boolean vector.
  #    - The boolean vector is result of single base pattern matching (e.g. G|C, G, etc.) 
  #      of a DNA sequence.
  #    - width is the size of the sliding windows.
  #    - Boolean vector is smaller in size compare to DNA sequence. e.g. chr1 size is 2 GB but
  #      it's G|C boolean vector is only 900 MB
  
  total.bins <- length(bool.vector) - width + 1
  
  if (ncpu == 1) {

    # initialise bool.percent vector
    bool.percent <- rep(NA, total.bins)

    # calculate the percentage
    i <- 1
    while ( i <= total.bins ) {
      bool.percent[i] <- sum(bool.vector[ i:(i+width-1) ]) / width * 100
      i <- i + 1
    }
    
  } else {
    
    # Distribute bins to ncpu
    chunks <- distributeChunk(total.bins, ncpu)

    bool.vector.chunk <- lapply(seq_along(chunks$start), function(i) bool.vector[ chunks$start[i]:(chunks$end[i]+width-1) ])

    # Calculate GC percentage
    bool.percent <- foreach(bool.vector=bool.vector.chunk, .combine='c', .noexport="bool.vector") %dopar% {
      
      total.bins <- length(bool.vector) - width + 1
      
      # initialise bool_percent vector
      bool.percent <- rep(NA, total.bins)
      
      # calculate the percentage
      i <- 1
      while ( i <= total.bins ) {
        bool.percent[i] <- sum(bool.vector[ i:(i+width-1) ]) / width * 100
        i <- i + 1
      }
      
      return(bool.percent)
    }
  }
  
  return(bool.percent)
}


# Binning method
# 1) Loop through 1:width
# 2) Form a matrix from a vector where column 1 is bin 1 and so forth.
# 3) Use colSums() to calculate percentage.
# 4) Sort so that it reflect a sliding window method.
#   -faster but not memory efficient - reached 12 GB for chr1 at 100 width
# 
## Initiate vector
#GC.percent <- rep(NA, length(bool.vector) - width + 1)
#
# for (i in 1:w) {
#   
#   # populate the GC.percent to indicate sliding. My method is binning from 1 to w
#   GC.percent[seq(i, length(GC.percent), width)] <- 
#     colSums( matrix( GC.bool[ i:( (chromosomes.size[[chr]]-i+1) %/% width * width + i -1) ], nrow = width) ) / width * 100
#   
# }