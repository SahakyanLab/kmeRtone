countSlidingBool <- function(bool.vector, width, ncpu=1) {
  # This function count the boolean "TRUE" of sliding windows of a boolean vector.
  #    - The boolean vector is result of single base pattern matching (e.g. G|C, G, etc.) 
  #      of a DNA sequence.
  #    - width is the size of the sliding windows.
  #    - Boolean vector is smaller in size compare to DNA sequence. e.g. chr1 size is 2 GB but
  #      it's G|C boolean vector is only 900 MB
  
  total.bins <- length(bool.vector) - width + 1
  
  if (ncpu == 1) {

    # initialise count.vector vector
    #count.vector <- rep(NA, total.bins)
    
    # initialise count table
    count.table <- rep(0, 101)
    names(count.table) <- 0:100

    # count the "TRUE"
    i <- 1
    while ( i <= total.bins ) {
      count.percent <- as.integer(sum(bool.vector[ i:(i+width-1) ]) / width * 100)
      count.table[count.percent+1] <- count.table[count.percent+1] + 1
      i <- i + 1
    }
    
  } else {
    
    # Distribute bins to ncpu
    chunks <- distributeChunk(total.bins, ncpu)

    bool.vector.chunk <- lapply(seq_along(chunks$start), function(i) bool.vector[ chunks$start[i]:(chunks$end[i]+width-1) ])

    toExport <- c("width")
    
    # Count the "TRUE"
    count.table <- foreach(bool.vector=bool.vector.chunk, .combine='c', .noexport=ls()[!ls() %in% toExport]) %dopar% {
      
      total.bins <- length(bool.vector) - width + 1
      
      # initialise count table
      count.table <- rep(0, 101)
      names(count.table) <- 0:100
      
      # count the "TRUE"
      i <- 1
      while ( i <= total.bins ) {
        count.percent <- as.integer(sum(bool.vector[ i:(i+width-1) ]) / width * 100)
        count.table[count.percent+1] <- count.table[count.percent+1] + 1
        i <- i + 1
      }
      return(count.table)
    }
  }
  
  count.table <- tapply(count.table, names(count.table), sum)
  count.table <- count.table[order(as.integer(names(count.table)))]
  
  return(count.table)
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