calKmerSkew <- function(kmer.table) {
  
  # kmer.table  <data.table>  Table with 3 columns: kmer, pos_strand, neg_strand
  
  kmer.table[, kmer_skew := (pos_strand - neg_strand) /
               (pos_strand + neg_strand)]
  
  
  invisible(kmer.table)
}