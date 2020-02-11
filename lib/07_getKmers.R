getKmers <- function(env){
  
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
  len <- env$genomic.coordinate[, unique(end - start + 1)]
  
  cat("--Length of genomic coordinates:", len, "nt\n")
  
  env$genomic.coordinate[, c("original_start", "original_end") := list(start, end)]
  
  # expand the coordinate to size k if it is smaller
  if (length(len) == 1 & len < env$k) {
    
    cat("--Length of the coordinates is smaller than size k.\n")
    cat("--Scaling up to size k.\n")
    
    # expand to k
    expandGenCoordinate("genomic.coordinate", env$kmertone.env, env$k)
    
    # trim down if out of range coordinate
    trimGenCoordinates("genomic.coordinate", "genome", remove = FALSE, env$kmertone.env)
    
    # flag region with less than size k (as a result of trimming down)
    env$genomic.coordinate[ end-start+1 < k, c("start", "end") := NA ]
    
  } else if (len >= env$k) {
    
    cat("--Coordinate lengths are the same or bigger than size k.\n--Kmer extraction will",
        "be performed by sliding window of size k.\n")
  }
  
  # remove overlapping case region
  cat("--Detecting overlapping case regions. ")
  before <- env$genomic.coordinate[(!is.na(start)) & (!is.na(end)), .N]
  
  removeAllOverlaps("genomic.coordinate", env, remove = FALSE)
  
  loss <- (before - env$genomic.coordinate[(!is.na(start)) & (!is.na(end)), .N]) / before * 100
  
  cat(loss, "% overlapping case regions are removed.\n")
  
  cat("Extracting case kmers...")
  # now that genomic coordinates are resolved, extract the case kmers
  case.kmers <- extractKmers("genomic.coordinate", "genome", env$k, env$DNA.pattern, env$kmertone.env)
  
  ## ------------------------------
  # INSENSITIVE STRAND MODE
  if (env$strand.mode == "insensitive") countReverseComplement("case.kmers")
  
  
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")
  
  # -------------------------------------------------------------------------------------------
  # SELECTION OF CASE REGIONS
  
  start.time <- Sys.time()
  
  cat("Selecting case zone...")
  
  # copy table
  case.region <- env$genomic.coordinate[, .(chromosome, start = original_start, end = original_end, strand)]
  
  ## increase to k if case region lower than k
  if (length(len) == 1 & len < env$k) {
    expandGenCoordinate("case.region", k=env$k)
  } else if ( length(len) > 1 && sum(len < env$k) != length(len) ) {
    stop("Some case regions have shorter size than k. I cannot decide on what to do.")
  }
  
  # add buffer region
  case.region[, `:=`(start = start - control.relative.position[1],
                     end = end + control.relative.position[1])]
  
  # trim out of range coordinates
  trimGenCoordinates("case.region", "genome", remove = TRUE)
  
  # merge overlapped regions
  mergeGenCoordinate("case.region")
  
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")
  
  # ------------------------------------------------------------------------------------------
  # SELECTION OF CONTROL REGIONS
  
  start.time <- Sys.time()
  
  cat("Selecting control regions...")
  
  # select control regions
  control.region <- env$genomic.coordinate[!is.na(end-start), .(
    start = c(start - control.relative.position[2], # [upstream] start
              end + control.relative.position[1]),  # (downstream) start
    end = c(start - control.relative.position[1],   # [upstream] end
            end + control.relative.position[2])     # (downstream) end
  ), by = .(chromosome, strand)]
  
  rm(genomic.coordinate, envir = env)
  gc()
  
  # trim out of range coordinates
  trimGenCoordinates("control.region", "genome", remove = TRUE)
  
  # merge overlapped regions
  mergeGenCoordinate("control.region")
  
  # remove overlapping case-zone portion (taking too long - need to optimise)
  removeCaseZone("control.region", "case.region", env$strand.mode)
  
  # clear up the memory because we just memory copy genomic coordinates
  gc()
  
  # remove region with lower than size k
  control.region <- control.region[end-start+1>=k]
  
  # Save the coordinates
  fwrite(control.region, "data/control_regions_coordinates.csv")
  
  # length
  control.region[, max(end-start+1)]
  
  # remove unneeded object
  rm(case.region)
  gc()
  
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")
  
  # ------------------------------------------------------------------------------------------
  # EXTRACTION OF CONTROL KMERS
  
  start.time <- Sys.time()
  
  cat("Extracting control kmers...\n")
  
  control.kmers <- extractKmers("control.region", "genome", env$k, env$DNA.pattern, env2 = env$kmertone.env)
  
  ## ------------------------------
  # INSENSITIVE STRAND MODE
  if (env$strand.mode == "insensitive") countReverseComplement("control.kmers")
  
  time.diff <- Sys.time() - start.time
  cat("Extracting control kmers...DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")
  
  # -----------------------------------------------------------------------------------------
  # KMERS TABLE
  
  # rename column count to either case or control
  setnames(control.kmers, "count", "control")
  setnames(case.kmers, "count", "case")
  
  # merge the table
  env$kmers <- merge(control.kmers, case.kmers, all = TRUE)
  
  # change NA to 0
  env$kmers[is.na(control), control := 0]
  env$kmers[is.na(case), case := 0]
  
  fwrite(kmers, "data/kmers.csv")
  
  # -----------------------------------------------------------------------------------------
  # MISC
  
  # restore back genomic coordinate table to original if changed
  if ("original_start" %in% colnames(env$genomic.coordinate)) {
    env$genomic.coordinate[, `:=`(start = original_start, end = original_end)]
    env$genomic.coordinate[, `:=`(original_start = NULL, original_end = NULL)]
  }
}