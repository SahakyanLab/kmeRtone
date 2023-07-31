#' Read FASTA files.
#'
#' @param fasta.file A path to a FASTA file.
#'
#' @return A sequence vector with header names
#'
#' @export
readFASTA <- function(fasta.file) {
  # Read fasta files and return sequence named with its header. Multiple headers
  # results in a vector of sequences

  fasta <- readLines(fasta.file)

  idx.head <- grep("^>", fasta)
  header <- fasta[idx.head]

  idx.start <- idx.head + 1
  idx.end <- c(idx.head[-1] - 1, length(fasta))
  seq <- sapply(seq_along(header), function(i) {
    paste(fasta[idx.start[i]:idx.end[i]], collapse = "")
  })

  names(seq) <- header

  return(seq)
}
