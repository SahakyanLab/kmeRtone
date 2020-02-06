filterTable <- function(env) {
  # remove sequence that don't match damage pattern
  # remove mitochondria
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, DNA.pattern, chromosome.names
  #     Packages          : data.table
  #     Functions         : -
  
  options( scipen=999)
  
  cat("\nSummary of data that are removed from the genomic coordinate table.\n")
  cat("-------------------------------------------------------------------\n")
  
  # damage pattern summary
  pattern.count <- data.table(env$genomic.coordinate[, table(sequence) ])
  
  # if there is other than pattern of interest
  if (nrow(pattern.count[sequence != env$DNA.pattern]) > 0) {
    fwrite(pattern.count, "damage_pattern_count.csv")
    # calculate percentage
    pattern.percent <- env$genomic.coordinate[sequence != env$DNA.pattern, table(sequence) 
                                          / nrow(env$genomic.coordinate) * 100]
    
    cat("\nThe following percentage of unwanted sequence pattern are removed:\n")
    print(pattern.percent)
    cat("\n")
  }
  
  cat("Different DNA pattern\t:", pattern.count[sequence != env$DNA.pattern, sum(N) / 
                                                  nrow(env$genomic.coordinate) * 100], "%\n")
  
  # mitochondria - c("chrM", "chrmt")
  mito <- c("chrM", "chrmt")
  mito.percent <- env$genomic.coordinate[chromosome %in% mito, .N] / nrow(env$genomic.coordinate) * 100
  
  cat("Mitochondria\t\t:", mito.percent, "%\n")
  
  genomic.coordinate.remove <- env$genomic.coordinate[chromosome %in% mito | sequence != env$DNA.pattern]
  fwrite(genomic.coordinate.remove, "data/removed_data_table.csv")
  
  cat("Total\t\t\t:", nrow(genomic.coordinate.remove)/nrow(env$genomic.coordinate)*100, "%\n")
  cat("-------------------------------------------------------------------\n\n")
  
  env$genomic.coordinate <- env$genomic.coordinate[(chromosome != "chrM") & (sequence == env$DNA.pattern)]
  env$chromosome.names <- env$chromosome.names[!env$chromosome.names %in% mito]

}