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
  
  # expand the coordinate to size k if it is smaller
  if (length(len) == 1 & len < env$k) {
    
    cat("--Length of the coordinates is smaller than size k.\n")
    cat("--Scaling up to size k.\n")
    
    # back up original coordinates
    env$genomic.coordinate[, `:=`(original_start = start, original_end = end)]
    
    # expand to k
    expandGenCoordinate("genomic.coordinate", env$kmertone.env, env$k)
    
    # trim down if out of range coordinate
    trimGenCoordinates("genomic.coordinate", "genome", trim = TRUE, env$kmertone.env)
    
  } else if (len >= env$k) {
    
    cat("--Coordinate lengths are the same or bigger than size k.\n--Kmer extraction will",
        "be performed by sliding window of size k.\n")
  }
  
  # now that genomic coordinates are resolved, extract the case kmers
  case.kmers <- extractKmers("genomic.coordinate", "genome", env$k, env$DNA.pattern, env$kmertone.env)
  
  ## ------------------------------
  # INSENSITIVE STRAND MODE
  if (env$strand.mode == "insensitive") countReverseComplement("case.kmers")
  
  
  time.diff <- Sys.time() - start.time
  cat("Extracting case kmers...DONE! ---", time.diff[1], attr(time.diff, "units"), "\n\n")
  
  # -------------------------------------------------------------------------------------------
  # SELECTION OF CASE REGIONS
  
  start.time <- Sys.time()
  
  cat("Selecting case zone...")
  
  # create case zone coordinate i.e. kmer/case region + buffer
  case.region <- env$genomic.coordinate[, .(chromosome, start, end, strand)]
  case.region[, `:=`(start = start - control.relative.position[1],
                   end = end + control.relative.position[1])]
  
  # trim out of range coordinates
  trimGenCoordinates("case.region", "genome", trim = TRUE)
  
  # merge overlapped regions
  mergeGenCoordinate("case.region")
  
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n\n")
  
  # ------------------------------------------------------------------------------------------
  # SELECTION OF CONTROL REGIONS
  
  start.time <- Sys.time()
  
  cat("Selecting control regions...")
  
  # select control regions
  control.region <- env$genomic.coordinate[, .(
    start = c(start - control.relative.position[2], # [upstream] start
              end + control.relative.position[1]),  # (downstream) start
    end = c(start - control.relative.position[1],   # [upstream] end
            end + control.relative.position[2])     # (downstream) end
  ), by = .(chromosome, strand)]
  
  # trim out of range coordinates
  trimGenCoordinates("control.region", "genome", trim = TRUE)
  
  # merge overlapped regions
  mergeGenCoordinate("control.region")
  
  # remove overlapping case-zone portion
  removeCaseZone("control.region", "case.region", env$strand.mode)
  
  # remove region with size less than 2
  control.region <- control.region[end - start + 1 > 2]
  
  # clear up the memory because we just memory copy genomic coordinates
  gc()
  
  # Save the coordinates
  fwrite(control.region, "data/control_regions_coordinates.csv")
  
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n\n")
  
  # ------------------------------------------------------------------------------------------
  # EXTRACTION OF CONTROL KMERS
  
  start.time <- Sys.time()
  
  cat("Extracting control kmers...\n")
  
  control.kmers <- extractKmers("control.region", "genome", env$k, env$DNA.pattern, env2 = env$kmertone.env)
  
  ## ------------------------------
  # INSENSITIVE STRAND MODE
  if (env$strand.mode == "insensitive") countReverseComplement("control.kmers")
  
  time.diff <- Sys.time() - start.time
  cat("Extracting control kmers...DONE! ---", time.diff[1], attr(time.diff, "units"), "\n\n")
  
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
  
  # -----------------------------------------------------------------------------------------
  # MISC
  
  # restore back genomic coordinate table to original if changed
  if ("original_start" %in% colnames(env$genomic.coordinate)) {
    env$genomic.coordinate[, `:=`(start = original_start, end = original_end)]
    env$genomic.coordinate[, `:=`(original_start = NULL, original_end = NULL)]
  }
}