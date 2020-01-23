checkCoordinate <- function(genomic.coordinate, genome, damage.pattern) {
  
  # Dependency: data.table
  
  dt <- fread(genomic.coordinate.path)
  
  # rename the column
  col.names = c("chromosome", "start", "end", "strand")
  if (colnames(dt[, 1:4]) != col.names) {
    for (i in seq_along(col.names)) {
      setnames(dt, colnames(dt)[i], col.names[i])
    }
  }
  
  # check that chromosome names in the table and genome are consistent
  if (sum( order(chromosome.names) %in% dt[, order(unique(chromosome))] ) == length(chromosome.names) ) {
    stop("Chromosome names are not consistent")
  }
  
  # check that the coordinates point to pattern
  checkPattern(genomic.coordinate, genome, damage.pattern)
  
  
  
  return(dt)
}