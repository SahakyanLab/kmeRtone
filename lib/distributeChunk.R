distributeChunk <- function(total, n) {
  # This function will distribute chunks from [total] equally [n] number.
  # The output would be start and end index of the chunks.
  # It is useful when to manually distribute work equally to CPU, hence reducing
  # the overhead and memory usage. Within foreach loop, there will be an internal
  # loop of the chunk.
  
  # Distribute total to n
  chunk <- rep(total %/% n, n)
  remainder <- total %% n
  if (remainder > 0) {
    remainder.parts <- rep(0, n)
    remainder.parts[1:remainder] <- rep(1, remainder)
    chunk <- chunk + remainder.parts
  }
  
  # if there is zero
  idx.zero <- chunk == 0
  
  start <- c(1, cumsum(chunk) + 1)[1:n]
  end <- cumsum(chunk)
  
  start <- start[!idx.zero]
  end <- end[!idx.zero]
  
  idx <- list(start=start, end=end, size=chunk)
  
  return(idx)
}