prepGenCoordinate <- function(env) {

  # Dependencies:
  #     Kmertone variables : genomic.coordinate, chromosome.names, strand.mode
  #     Package            : data.table
  
  # processing a string format
  if (class(env$genomic.coordinate)[1] == "character") {

    filename <- env$genomic.coordinate
    
    # genomic.coordinate can be replicates
    env$genomic.coordinate <- lapply(env$genomic.coordinate, fread)
    
    if (sum(grepl("\\.bed", filename)) == length(filename)) {
      
      for (i in seq_along(env$genomic.coordinate)) {
        
        if (length(colnames(env$genomic.coordinate[[i]])) >= 6) {
          # strand sensitive i.e. has strand information
          
          env$genomic.coordinate[[i]][, colnames(env$genomic.coordinate[[i]])[!colnames(env$genomic.coordinate[[i]]) %in% 
                                                                                c("V1", "V2", "V3", "V6")] := NULL]
          
          setnames(env$genomic.coordinate[[i]], c("V1", "V2", "V3", "V6"), c("chromosome", "start", "end", "strand"))
          
          # bed use zero-based index and open-end index e.g. [0,10) which means from index 0 until index 9
          # change to R indexing
          env$genomic.coordinate[[i]][, start := start + 1]
          
        } else {
          
          if (strand.mode == "sensitive") {
            stop("You specify strand sensitive but there is no strand information in the table.")
          } else {
            
            # strand insensitive i.e. has no strand information
            env$genomic.coordinate[[i]][, colnames(env$genomic.coordinate[[i]])[!colnames(env$genomic.coordinate[[i]]) %in% 
                                                                                  c("V1", "V2", "V3")] := NULL]
            
            setnames(env$genomic.coordinate[[i]], c("V1", "V2", "V3"), c("chromosome", "start", "end"))
            env$genomic.coordinate[[i]][, strand := "*"]
          }
        }
      }
    }
  }
  
  # add column replicate for list of genomic.coordinate replicates
  if (class(env$genomic.coordinate)[1] == "list") {
    
    cat("Detecting", length(env$genomic.coordinate), "replicates\n")
    
    # rename columns
    for (rep in env$genomic.coordinate) {
      setnames(rep, colnames(rep)[1:4], c("chromosome", "start", "end", "strand"))
    }
    
    # add column replicate_number
    cnt <- 1
    for (i in 1:length(env$genomic.coordinate)) {
      env$genomic.coordinate[[i]][, paste0("replicate_", cnt) := 1]
      cnt <- cnt + 1
    }

    cat("Merging the genomic coordinate tables...")
    dt <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = c("chromosome", "start", "end", "strand"),
                                          all = TRUE),
                 x = env$genomic.coordinate)
    
    dt <- unique(dt)
    
    # reassign the consolidated table to genomic.coordinate
    env$genomic.coordinate <- dt
    
    # memory copy, so gc() to clear memory
    gc()
    
    cat("DONE\n")
    
  } else {
    
    setnames(x = env$genomic.coordinate, 
             old = colnames(env$genomic.coordinate)[1:4],
             new = c("chromosome", "start", "end", "strand"))
    
  }
  
  setkey(env$genomic.coordinate, chromosome)
  if (sum(!env$genomic.coordinate[, unique(chromosome)] %in% env$chromosome.names) > 0) {
    
    message(paste0("Genome chromosome names: ", paste(env$chromosome.names, collapse = " ")))
    message(paste0("Genomic coordinate table chromosome names: ", paste(env$genomic.coordinate[, unique(chromosome)],
                                                                        collapse = " ")))
    stop("Chromosome names are not consistent between genome and genomic coordinate table")
    
  }
}
