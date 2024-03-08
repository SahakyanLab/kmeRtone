#' Function calculates scores for k-mers based on case and control k-mer counts.
#' 
#' @param case.kmers A data.table containing k-mer counts in case samples.
#' @param control.kmers A data.table containing k-mer counts in control samples.
#' 
#' @return A data.table containing scores for each k-mer.
#' 
#' @importFrom data.table set setcolorder setkeyv setnames
#' 
#' @export
getScores <- function(case.kmers, control.kmers) {
  
  # Calculate k-mer skew
  #if (!strand.sensitive) {
    calKmerSkew(case.kmers)
    calKmerSkew(control.kmers)
  #}
  
  # Sum counts from both strands
  case.kmers[, case := pos_strand + neg_strand]
  control.kmers[, control := pos_strand + neg_strand]
  
  # Remove separate strand counts
  case.kmers[, c("pos_strand", "neg_strand") := NULL]
  control.kmers[, c("pos_strand", "neg_strand") := NULL]
  
  # Rename kmer_skew
  setnames(case.kmers, "kmer_skew", "case_skew")
  setnames(control.kmers, "kmer_skew", "control_skew")
  
  setkeyv(case.kmers, names(case.kmers)[grepl("kmer", names(case.kmers))])
  setkeyv(case.kmers, names(control.kmers)[grepl("kmer", names(case.kmers))])
  
  kmer.table <- merge(case.kmers, control.kmers, all = TRUE)
  
  scoreKmers(kmer.table)
  
  setcolorder(kmer.table,
              c(names(kmer.table)[grepl("kmer", names(kmer.table))],
                "case", "control", "case_skew", "control_skew", "z"))
  
  return(kmer.table)
}