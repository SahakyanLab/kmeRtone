#' Function aims to determine the optimum k-mer size for genomic analysis
#' based on the provided case data and genome. It may involve analyzing central
#' patterns within the genomic sequences.
#'
#' @param case Coordinate class object or similar structure for case data.
#' @param genome Genome class object or similar structure.
#' @param central.pattern Central pattern of the k-mer.
#'
#' @return Optimum k-mer size for the given genomic data.
GET_OPTIMUM_K <- function(case, genome, central.pattern) {
  
  # options(warn = 2); 
  library(kmeRtone)
  genome <- loadGenome("hg19", "UCSC")
  case <- loadCoordinate("../../../09_Hilary/github/kmer_scores/data/damage/CT_CPD/hg19_1/coordinate/csv/", single.len = 2, merge.replicates = T, rm.dup = T)
  
  for (chr.name in genome$chr_names) {
    
  }
}