filterTable <- function() {
  # remove sequence that don't match damage pattern
  # remove mitochondria
  
  # Dependencies:
  #     Global variables: genomic.coordinate, DNA.pattern
  #     Packages        : data.table
  #     Functions       : -
  
  options( scipen=999)
  
  cat("\nSummary of data that are removed from the genomic coordinate table.\n")
  cat("-------------------------------------------------------------------\n")
  
  # damage pattern summary
  pattern.count <- data.table(genomic.coordinate[, table(sequence) ])
  
  # if there is other than pattern of interest
  if (nrow(pattern.count[sequence != DNA.pattern]) > 0) {
    fwrite(pattern.count, "damage_pattern_count.csv")
    # calculate percentage
    pattern.percent <- genomic.coordinate[sequence != DNA.pattern, table(sequence) / nrow(genomic.coordinate) * 100]
    
    cat("\nThe following unwanted sequence pattern are removed:\n")
    print(pattern.percent)
    cat("\n")
  }
  
  cat("Different DNA pattern\t:", pattern.count[sequence != DNA.pattern, sum(N) / nrow(pattern.count) * 100], "%\n")
  
  # mitochondria - c("chrM", "chrmt")
  mito <- c("chrM", "chrmt")
  mito.percent <- genomic.coordinate[chromosome %in% mito, .N] / nrow(genomic.coordinate) * 100
  
  cat("Mitochondria\t\t:", mito.percent, "%\n")
  
  genomic.coordinate.remove <- genomic.coordinate[chromosome %in% mito | sequence != DNA.pattern]
  fwrite(genomic.coordinate.remove, "data/removed_data_table.csv")
  
  cat("Total\t\t\t:", nrow(genomic.coordinate.remove)/nrow(genomic.coordinate)*100, "%\n")
  cat("-------------------------------------------------------------------\n")
  
  genomic.coordinate <<- genomic.coordinate[(chromosome != "chrM") & (sequence == DNA.pattern)]
  chromosome.names <<- chromosome.names[!chromosome.names %in% mito]

}