getControlKmers <- function(genomic.coordinate, genome, k, DNA.pattern, strand.mode, 
                            case.kmers, env=parent.frame()) {
  # genomic.coordinate needs to go through getCaseKmers function first where
  # overlapping kmers coordinates are flagged with NA and there are two additional
  # columns: original_start and original_end
  
  # Dependencies
  #    Kmertone variables  : genomic.coordinate 
  #    Functions           : extractKmers, reverseComplement
  
  # -------------------------------------------------------------------------------------------
  # SELECTION OF CASE REGIONS
  
  start.time <- Sys.time()
  
  cat("Selecting case zone...")
  
  # copy table
  case.region <- env[[genomic.coordinate]][, .(chromosome, start = original_start, end = original_end, strand)]
  
  ## increase to k if case region lower than k
  if (length(len) == 1 & len < k) {
    expandGenCoordinate("case.region", k=k)
  } else if ( length(len) > 1 && sum(len < k) != length(len) ) {
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
  control.region <- env[[genomic.coordinate]][!is.na(end-start), .(
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
  removeCaseZone("control.region", "case.region", strand.mode)
  
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
  
  control.kmers <- extractKmers("control.region", "genome", k, DNA.pattern, env2 = env)
  
  ## ------------------------------
  # INSENSITIVE STRAND MODE
  if (strand.mode == "insensitive") countReverseComplement("control.kmers")
  
  time.diff <- Sys.time() - start.time
  cat("Extracting control kmers...DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")
  
  # -----------------------------------------------------------------------------------------
  # KMERS TABLE
  
  # rename column count to either case or control
  setnames(control.kmers, "count", "control")
  setnames(env[[case.kmers]], "count", "case")
  
  # merge the table
  kmers <- merge(control.kmers, env[[case.kmers]], all = TRUE)
  
  # change NA to 0
  kmers[is.na(control), control := 0]
  kmers[is.na(case), case := 0]
  
  # -----------------------------------------------------------------------------------------
  # MISC
  
  # restore back genomic coordinate table to original if changed
  if ("original_start" %in% colnames(env[[genomic.coordinate]])) {
    env[[genomic.coordinate]][, `:=`(start = original_start, end = original_end)]
    env[[genomic.coordinate]][, `:=`(original_start = NULL, original_end = NULL)]
  }
  
  return(kmers)
}