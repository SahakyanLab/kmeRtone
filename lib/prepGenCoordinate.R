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
  
  # Processing a data.table format
  if (env$strand.mode == "sensitive") {
    col.names <- c("chromosome", "start", "end", "strand")
  } else if (env$strand.mode == "insensitive") {
    col.names <- c("chromosome", "start", "end")
  }
  
  # add column replicate for list of genomic.coordinate replicates
  if (class(env$genomic.coordinate)[1] == "list") {
    
    cat("Detecting", length(env$genomic.coordinate), "replicates\n")
    
    # rename columns
    for (rep in env$genomic.coordinate) {
      setnames(rep, colnames(rep), col.names)
    }
    
    # add column replicate_number
    cnt <- 1
    for (i in 1:length(env$genomic.coordinate)) {
      env$genomic.coordinate[[i]][, paste0("replicate_", cnt) := 1]
      cnt <- cnt + 1
    }

    cat("Merging the genomic coordinate tables...")
    dt <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = col.names, all = TRUE),
                 x = env$genomic.coordinate)
    
    dt <- unique(dt)
    
    # reassign the consolidated table to genomic.coordinate
    env$genomic.coordinate <- dt
    
    # memory copy, so gc() to clear memory
    gc()
    
    cat("DONE\n")
    
  } else {
    
    setnames(x = env$genomic.coordinate, 
             old = colnames(env$genomic.coordinate),
             new = col.names)
    
  }
  
  setkey(env$genomic.coordinate, chromosome)
  if (sum(!env$genomic.coordinate[, unique(chromosome)] %in% env$chromosome.names) > 0) {
    
    message(paste0("Genome chromosome names: ", env$chromosome.names))
    message(paste0("Genomic coordinate table chromosome names: ", env$genomic.coordinate[, unique(chromosome)]))
    stop("Chromosome names are not consistent between genome and genomic coordinate table")
    
  }
}
