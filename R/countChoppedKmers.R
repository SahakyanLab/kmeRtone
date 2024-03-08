#' Function chops k-mers within specified ranges of a sequence and counts 
#' them. It uses either a substring method or functionalities from the 
#' Biostrings package.
#' 
#' @param dna.seq A string of sequence.
#' @param starts Start positions.
#' @param ends End positions.
#' @param k Size of kmer.
#' @param method Method: "Biostrings" or "substring". Default is Biostrings.
#' @return A k-mer-named vector of counts.
#'
#' @importFrom stringi stri_sub
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' 
#' @export
countChoppedKmers <- function(dna.seq, starts, ends, k, method="auto") {
  
  if (length(dna.seq) > 1) stop("Please input only a single string.")
  
  if (method == "substring" | k > 12) {
    
    froms <- lapply(which(ends - starts + 1 >= k), function(i) {
      starts[i]:(ends[i] - k + 1)
    }) |> unlist(use.names = FALSE)
    
    counts <- stri_sub(dna.seq, froms, length = k) |> table() |> c()
    
    # Reformat the return to be the same like Biostrings if no k-mer found.
    if (length(counts) == 0) names(counts) <- character()

  } else if (method == "Biostrings" | k <= 12) {
    
    counts <- stri_sub(dna.seq, starts, ends) |> Biostrings::DNAStringSet() |>
      Biostrings::oligonucleotideFrequency(width = k, simplify.as = "collapsed")
    counts <- counts[counts > 0]
    
  }

  return(counts)
}