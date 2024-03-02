GET_OPTIMUM_K <- function(case, genome, central.pattern) {
  
  options(warn = 2); library(kmeRtone)
  genome <- loadGenome("hg19", "UCSC")
  case <- loadCoordinate("../../../09_Hilary/github/kmer_scores/data/damage/CT_CPD/hg19_1/coordinate/csv/", single.len = 2, merge.replicates = T, rm.dup = T)
  
  for (chr.name in genome$chr_names) {
    
  }
  
  # Plot
  
}