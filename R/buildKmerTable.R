#' Count k-mers from given sequence(s) and build a data.table of k-mer counts.
#'
#' Only existed k-mers are returned in data.table object.
#'
#' @param dna.seqs String of sequence(s).
#' @param k Size of kmer.
#' @param remove.N Remove unknown base? Default is TRUE.
#' @param method K-mer counting method: "Biostrings", "sliding", or "auto".
#'    Default is "auto"; For k > 8, sliding method is used.
#' @return A `data.table` object with column kmer and N.
#'
#' @export
buildKmerTable <- function(dna.seqs, k, method="auto", remove.N=TRUE) {
  
  if (!remove.N & method == "Biostrings") {
    stop("Chopping method Biostrings remove N automatically.")
  }
  
  if (method == "Biostrings" | (k <= 12 & remove.N)) {
    
    counts <- Biostrings::oligonucleotideFrequency(
      Biostrings::DNAStringSet(unlist(dna.seqs)),
      width = k,
      simplify.as = "collapsed")
    counts <- data.table(kmer = names(counts), N = counts)[N > 0]
    
  } else if (method %in% c("sliding", "auto")) {
    
    counts <- countKmers(unlist(dna.seqs), k)
    counts <- data.table(kmer = names(counts), N = counts)
    
    if (remove.N) counts <- counts[!kmer %like% "[^ACGT]"]
    
    setorder(counts, kmer)
    
  }
  
  return(counts)
}