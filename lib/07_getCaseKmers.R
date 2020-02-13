getCaseKmers <- function(genomic.coordinate, genome, k, DNA.pattern, strand.mode, 
                         remove.overlaps=TRUE, env=parent.frame()){
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, DNA.pattern, k, strand.mode
  #                         kmertone.env, genome
  #     Packages          : data.table, stringi
  #     Functions         : reverseComplement
  
  # Output: kmertone.env$kmers
  
  # -------------------------------------------------------------------------------------------
  # EXTRACTION OF CASE KMERS 
  
  start.time <- Sys.time()
  cat("Extracting case kmers...\n")
  
  # check if genomic coordinate length varies
  len <- env[[genomic.coordinate]][, unique(end - start + 1)]
  
  cat("--Length of genomic coordinates:", len, "nt\n")
  
  # expand the coordinate to size k if it is smaller
  if (length(len) == 1 & len < k) {
    
    cat("--Length of the coordinates is smaller than size k.\n")
    cat("--Scaling up to size k.\n")
    
    # expand to k
    expandGenCoordinate("genomic.coordinate", env, k)
    
    # trim down if out of range coordinate
    trimGenCoordinates("genomic.coordinate", genome, remove = FALSE, env)
    
    # flag region with less than size k (as a result of trimming down)
    env[[genomic.coordinate]][ end-start+1 < k, c("start", "end") := NA ]
    
  } else if (len >= k) {
    
    cat("--Coordinate lengths are the same or bigger than size k.\n--Kmer extraction will",
        "be performed by sliding window of size k.\n")
  }
  
  # remove overlapping case region
  if (remove.overlaps) {
    
    cat("--Detecting overlapping case regions. ")
    before <- env[[genomic.coordinate]][(!is.na(start)) & (!is.na(end)), .N]
    
    removeAllOverlaps("genomic.coordinate", env, remove = FALSE)
    
    loss <- (before - env[[genomic.coordinate]][(!is.na(start)) & (!is.na(end)), .N]) / before * 100
    
    cat(loss, "% overlapping case regions are removed.\n")
    
  } else {
    cat("Overlapping case kmers are set not to be removed.\n")
  }
  
  
  cat("Extracting case kmers...")
  # now that genomic coordinates are resolved, extract the case kmers
  case.kmers <- extractKmers("genomic.coordinate", genome, k, DNA.pattern, env)
  
  ## ------------------------------
  # INSENSITIVE STRAND MODE
  if (strand.mode == "insensitive") countReverseComplement("case.kmers")
  
  
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")
  
  return(case.kmers)
}