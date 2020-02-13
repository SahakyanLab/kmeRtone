filterTable <- function(genomic.coordinate, DNA.pattern, strand.mode, env=parent.frame()) {
  # remove sequence that don't match case pattern
  # remove mitochondria
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, DNA.pattern
  #     Packages          : data.table
  #     Functions         : -
  
  options( scipen=999)
  
  cat("\n\nSummary of data that are removed from the genomic coordinate table.\n")
  cat("-------------------------------------------------------------------\n")
  
  # ---------------------------------------------------------------------------------------
  # Different DNA pattern.
  if (!is.null(DNA.pattern)) {
    
    # case pattern summary
    pattern.count <- data.table(env[[genomic.coordinate]][, table(sequence) ])
    #fwrite(pattern.count, "data/case_pattern_count.csv")
    
    if (nrow(pattern.count[(!sequence %in% DNA.pattern)]) > 0) {
      
      # calculate percentage
      pattern.percent <- env[[genomic.coordinate]][(!sequence %in% DNA.pattern), table(sequence) 
                                                / nrow(env[[genomic.coordinate]]) * 100]
      
      cat("\nThe following percentage of unwanted sequence pattern are removed:\n")
      print(pattern.percent)
      cat("\n")
      
      diff.DNA.percent <- pattern.count[(!sequence %in% DNA.pattern), sum(N) / 
                                          nrow(env[[genomic.coordinate]]) * 100]
      
    } else {
      diff.DNA.percent <- 0
    }
    
  } else if (is.null(DNA.pattern)) {
    diff.DNA.percent <- NA
  }
  
  # ---------------------------------------------------------------------------------------
  # Mitochondria - c("chrM", "chrmt")
  mito <- c("chrM", "chrmt", "chrMT", "chrmito", "chrMito")
  if (env[[genomic.coordinate]][chromosome %in% mito, .N] > 0) {
    mito.percent <- env[[genomic.coordinate]][chromosome %in% mito, .N] / nrow(env[[genomic.coordinate]]) * 100
  } else {
    mito.percent <- 0
  }
  
  # ---------------------------------------------------------------------------------------    
  # If case happen on both strand in sensitive mode
  if (strand.mode == "sensitive" & (env[[genomic.coordinate]][strand == "*", .N] > 0 | 
                                    sum(duplicated(env[[genomic.coordinate]][, .(chromosome, start, end)])) > 0) ) {
    message("The strand mode is sensitive but there is information on both strand.")
    
    dup.idx <- duplicated(env[[genomic.coordinate]][, .(chromosome, start, end)]) |
      duplicated(env[[genomic.coordinate]][, .(chromosome, start, end)], fromLast = TRUE)
    
    dup.dt <- env[[genomic.coordinate]][dup.idx]
    setkey(dup.dt, chromosome, start, end)
    
    print(dup.dt)
    cat("\nBelow are the percentage of such incidence.\n")
    print(dup.dt[, table(sequence) / nrow(env[[genomic.coordinate]]) * 100])
    cat("\n")
    
    both.strand.percent <- dup.dt[, .N / nrow(env[[genomic.coordinate]]) * 100]
    
  } else if (strand.mode == "insensitive") {
    both.strand.percent <- NA
    dup.idx <- FALSE
  } else {
    both.strand.percent <- 0
    dup.idx <- FALSE
  }
  
  # Save removed data
  if (!is.null(DNA.pattern)) {
    total.remove <- env[[genomic.coordinate]][dup.idx | chromosome %in% mito | (!sequence %in% DNA.pattern), .N]
    #fwrite(genomic.coordinate[dup.idx | chromosome %in% mito | (!sequence %in% env$DNA.pattern)],
    #       "data/removed_data_table.csv")
  } else if (is.null(DNA.pattern)) {
    total.remove <- env[[genomic.coordinate]][dup.idx | chromosome %in% mito, .N]
    #fwrite(genomic.coordinate[dup.idx | chromosome %in% mito],
    #       "data/removed_data_table.csv")
  }
  
  total <- nrow(env[[genomic.coordinate]])
  # Remove and update the table
  if (!is.null(DNA.pattern)) {
    env[[genomic.coordinate]] <- env[[genomic.coordinate]][(!dup.idx) & (!chromosome %in% mito) & (sequence %in% DNA.pattern)]
  } else if (is.null(DNA.pattern)) {
    env[[genomic.coordinate]] <- env[[genomic.coordinate]][(!dup.idx) & (!chromosome %in% mito)]
  }
  gc()
  
  cat("Different DNA pattern\t:", diff.DNA.percent, "%\n")
  cat("Mitochondria\t\t:", mito.percent, "%\n")
  cat("Case on both strand\t:", both.strand.percent, "%\n")
  cat("Total\t\t\t:", total.remove/total*100, "%\n")
  cat("-------------------------------------------------------------------\n\n")
}