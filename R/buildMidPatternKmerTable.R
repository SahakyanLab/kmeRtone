#' Count k-mers with specified middle pattern from given sequence(s) and build
#'    a data.table of k-mer counts.
#'
#' Only existed k-mers are returned in data.table object.
#'
#' @param dna.seqs String of sequence(s).
#' @param k Size of kmer.
#' @param mid.patterns Middle patterns.
#' @param remove.N Remove unknown base? Default is TRUE.
#' @return A `data.table` object with column kmer and N.
#'
#' @export
buildMidPatternKmerTable <- function(dna.seqs, k, mid.patterns, remove.N=TRUE) {
  
  counts <- lapply(mid.patterns, function(mid.pattern)
    countMidPatternKmers(unlist(dna.seqs), k, mid.pattern)) |>
    unlist(recursive = FALSE)
  
  counts <- data.table(kmer = names(counts), N = counts)
    
  if (remove.N) counts <- counts[!kmer %like% "[^ACGT]"]
  
  setorder(counts, kmer)
  
  return(counts)
}