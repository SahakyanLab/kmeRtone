checkCoordinate <- function(genomic.coordinate, genome.path, damage.pattern, ncpu=1) {
  
  # Dependency: data.table
  
  dt <- genomic.coordinate
  
  # ------------ column names -----------------------------------------------------------
  # column names must be chromosome, start, end, and strand
  
  col.names = c("chromosome", "start", "end", "strand")
  
  if ( sum(!colnames(dt) %in% col.names) > 0 ) {
    stop("Please rename your columns to chromosome, start, end, strand")
  }
  
  # -------- chromosome names -----------------------------------------------------------
  # check that chromosome names in the table and genome are consistent
  
  chromosomes.filename <- list.files(genome.path)
  chromosome.names <- sub("\\..*", "", chromosomes.filename)
  
  chromosome.in.table <- dt[, unique(chromosome)]
  
  if ( sum(!chromosome.in.table %in% chromosome.names ) > 0 ) {
    message("Chromosome in the table:")
    print(chromosome.in.table)
    message("Chromosome in the genome:")
    print(chromosome.names)
    message("Please rename the chromosome appropriately.")
    stop("Chromosome names in genome and genomic coordinate table are not consistent!")
  }
  
  # ---------- damage coordinate ---------------------------------------------------------
  # check that the coordinates point to pattern
  
  # add column sequence to the table
  addColumnSequence(dt, genome.path, ncpu = 12)

  damage.summary <- dt[, table(sequence)]
  if (damage.summary[[damage.pattern]] != 100) {
    message("There is other pattern in the table:")
    print(damage.summary)
  }
  
  checkPattern(genomic.coordinate, genome, damage.pattern)
  
  return(dt)
}