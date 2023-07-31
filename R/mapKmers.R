#' Map k-mers of a given sequence and coordinate.
#' 
#' @field seq A single sequence string.
#' @field start A single of vector of start coordinates where k-mers are mapped.
#' @field end A single of vector of end coordinates where k-mers are mapped. It
#'      is NULL when all regions have the same length.
#' @field len A single integer for specifying fixed length of regions. The end
#'      is NULL in this case.
#' @field k Length of k-mer to be mapped.
#' @field rm.trunc.kmer Remove truncated k-mers as a result of out-of-bound
#'      region. Default is TRUE.
#'
#' @return Mapped k-mers.
#'
#' @export
mapKmers <- function(seq, start, end=NULL, len=NULL, k, rm.trunc.kmer=TRUE) {
  
  # Map sequence based on these following possible situations:
  
  # Situation:
  # 1. Only start position: Exact k-mer extract (coordinate in k-mer state)
  # 2. Only start and end exist i.e. varied region length: sliding window
  # 3. Only start and len exist i.e. big fixed region (len > k): sliding window
  
  # Note:
  # This function primarily relies on stri_sub function which has a special
  #     meaning for negative position i.e. starting from the back. So this
  #     function avoids that by removing or trimming negative or zero value
  #     to 1.
  
  # Error checking
  if(!is.null(end) & !is.null(len)) {
    stop("Either end or len should be input.")
  }
  
  do.sliding.windows <- !is.null(end) | (!is.null(len) && len > k)
  do.exact.mapping <- is.null(end) & is.null(len)
  
  # Remove region less than k
  if (!is.null(end)) {
    idx <- end - start + 1 >= k
    start <- start[idx]
    end <- end[idx]
  }
  
  # Resolve k-mer start positions for sliding big region
  if (do.sliding.windows) {
    
    start <- unlist(
      Map(function(start, end) start:(end - k + 1),
          start,
          if(!is.null(end)) end else start + len - 1))
    
  }
  
  # Resolve truncated k-mers because of negative or zero value
  if (rm.trunc.kmer) {

    start <- start[start >= 1]
    
  } else if (!rm.trunc.kmer) {
    
    start[(start + k - 1) < 1] <- NA
    idx.trunc <- start < 1 & (start + k - 1) > 1
    idx.trunc[is.na(idx.trunc)] <- FALSE
    start.trunc <- start[idx.trunc]
    start[idx.trunc] <- 1
    
  }
  
  # Mapping operation
  kmers <- stri_sub(seq, start, length = k)
  
  # Remove truncated k-mers as a result of out-of-bound far right region.
  if (rm.trunc.kmer) {
    
    kmers <- kmers[stri_length(kmers) == k]
    
  } else if (!rm.trunc.kmer) {
    
    kmers[idx.trunc] <- stri_sub(kmers[idx.trunc],
                                 from = 1,
                                 length = k + start.trunc - 1)
    
  }
  
  return(kmers)
}