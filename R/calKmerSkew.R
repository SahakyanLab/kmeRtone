#' Function calculates the skew of k-mers based on their occurrence in 
#' positive and negative strands.
#'
#' @param kmer.table data.table with columns: kmer, pos_strand, neg_strand.
#' @return data.table with the kmer_skew column.
#'
#' @export
calKmerSkew <- function(kmer.table) {
  
  # kmer.table  <data.table>  Table with 3 columns: kmer, pos_strand, neg_strand
  
  kmer.table[, kmer_skew := (pos_strand - neg_strand) /
               (pos_strand + neg_strand)]
  
  
  invisible(kmer.table)
}