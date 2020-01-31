countWindowBoolPercentage <- function(bool.vector, start, end, ncpu=1) {
  # This function calculate percentage of TRUE of a given bool.vector window.
  # 
  # bool.vector    a vector of boolean e.g. c(TRUE, FALSE, FALSE)
  # start          start position of a window. It can be a vector
  # end            end position of a window. It can a vector
  
  # Dependencies: distributeChunk(), parallel, doParallel
  
  if (ncpu == 1) {
    
    count.vector <- mapply(function(start, end) sum(bool.vector[start:end]) / (end-start+1) * 100,
                           start, end)
    
  } else {
    
    chunks <- distributeChunk(length(start), ncpu)
    
    start.chunk <- lapply(1:ncpu, function(i) start[ chunks$start[i]:chunks$end[i] ])
    end.chunk <- lapply(1:ncpu, function(i) end[ chunks$start[i]:chunks$end[i] ])
    
    bool.vector.chunk <- lapply(1:ncpu, function(i) bool.vector[ min(start.chunk[[i]]):max(end.chunk[[i]]) ])
    
    count.vector <- foreach(start=start.chunk, end=end.chunk, bool.vector=bool.vector.chunk, 
                            .noexport = "bool.vector", .combine = "c") %dopar% {
      
      # offset the coordinate
      start.min <- min(start)
      start <- start - start.min + 1
      end <- end - start.min + 1
      
      count.vector <- mapply(function(start, end) sum(bool.vector[start:end]) / (end-start+1) * 100,
                             start, end)
      
      return(count.vector)
    }
  }
  
  return(count.vector)
}