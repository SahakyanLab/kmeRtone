prepGenCoordinate <- function() {
  # global variable used: genomic.coordinate, chromosome.names

  # add column replicate for list of genomic.coordinate replicates
  if (class(genomic.coordinate)[1] == "list") {
    
    cat("Detecting", length(genomic.coordinate), "replicates\n")
    
    # rename columns
    for (rep in genomic.coordinate) {
      setnames(rep, colnames(rep)[1:4], c("chromosome", "start", "end", "strand"))
    }
    
    # add column replicate_number
    cnt <- 1
    for (i in 1:length(genomic.coordinate)) {
      genomic.coordinate[[i]][, paste0("replicate_", cnt) := 1]
      cnt <- cnt + 1
    }

    cat("Merging the genomic coordinate tables...")
    dt <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = c("chromosome", "start", "end", "strand"), all = TRUE),
                 x = genomic.coordinate)
    
    dt <- unique(dt)
    
    # reassign the consolidated table to genomic.coordinate
    genomic.coordinate <<- dt
    
    cat("DONE\n")
    
  } else {
    
    setnames(genomic.coordinate, colnames(genomic.coordinate), c("chromosome", "start", "end", "strand"))
    
  }
  
  setkey(genomic.coordinate, chromosome)
  if (sum(!genomic.coordinate[, unique(chromosome)] %in% chromosome.names) > 0) {
    
    error("Chromosome names are not consistent between genome and genomic coordinate table")
    message(paste0("Genome chromosome names: ", chromosome.names))
    message(paste0("Genomic coordinate table chromosome names: ", genomic.coordinate[, unique(chromosome)]))
    
  }
  
}