filterTable <- function(genomic.coordinate, damage.pattern) {
  # remove sequence that don't match damage pattern
  # remove mitochondria
  
  dt <- genomic.coordinate
  
  # damage pattern summary
  pattern.count <- data.table(dt[, table(sequence) ])
  
  # calculate percentage
  pattern.count[, N / sum(N) * 100]
  
  fwrite(pattern.count, "damage_pattern_count.csv")
  
  
  # mitochondria
  dt[chromosome == "chrM", .N] / nrow(dt) * 100
  
  dt.remove <- dt[chromosome == "chrM"]
  dt.remove <- rbind(dt.remove, dt[sequence != damage.pattern])
  
  fwrite(dt.remove, "data/filter/removed_data.csv")
  
  dt <- dt[chromosome != "chrM" | sequence != damage.pattern]
  
  return(dt)
}