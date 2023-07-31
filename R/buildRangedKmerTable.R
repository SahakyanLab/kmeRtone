#' Count kmers from a sequence in given ranges and build a data.table of k-mer
#'    counts.
#' 
#' @param dna.seq String of sequence.
#' @param starts Start positions.
#' @param ends End positions.
#' @param k Size of kmer.
#' @param method Method options: "sliding" or "chopping". Chopping consumes a
#'    lot of memory for extremely long sequence using "substring" method. Using
#'    "Biostrings" for k > 12 is memory consuming. Default is "sliding".
#' @param chopping.method Chopping method: "Biostrings" or "substring". Default
#'    is "auto". 
#' @param remove.N Remove unknown base N? Default is TRUE.
#' @return A `data.table` object with column kmer and N.
#'
#' @export
buildRangedKmerTable <- function(dna.seq, starts, ends, k, method="sliding",
                                chopping.method="auto", remove.N=TRUE) {
  
  if (!remove.N & chopping.method == "Biostrings") {
    stop("Chopping method Biostrings remove N automatically.")
  }
  
  if (method == "sliding") {
    counts <- countRangedKmers(unlist(dna.seq), starts, ends, k)
  } else if (method == "chopping") {
    if (!remove.N & chopping.method == "auto") chopping.method <- "substring"
    counts <- countChoppedKmers(dna.seq, starts, ends, k, chopping.method)
  }
  
  counts <- data.table(kmer = names(counts), N = counts)
  
  if (remove.N) counts <- counts[!kmer %like% "[^ACGT]"]
  
  setorder(counts, kmer)
  
  return(counts)
}