getControlKmers <- function(genomic.coordinate, genome,
                            control.relative.position, k, DNA.pattern,
                            strand.sensitive, case.kmers,
                            env = parent.frame()) {
  # genomic.coordinate needs to go through getCaseKmers function first where
  # overlapping kmers coordinates are flagged with NA and there are two
  # additional columns: original_start and original_end

  # Dependencies
  #    Kmertone variables  : genomic.coordinate
  #    Functions           : extractKmers, reverseComplement

  # ----------------------------------------------------------------------------
  # SELECTION OF CASE REGIONS

  start.time <- Sys.time()

  cat("Selecting case zone...")

  # copy table
  case.region <- env[[genomic.coordinate]][, .(chromosome, start, end, strand)]

  # # check if genomic coordinate length varies
  # len <- case.region[, unique(end - start + 1)]
  #
  # ## increase to k if case region lower than k
  # if (length(len) == 1 & len < k) {
  #   expandGenCoordinate("case.region", k=k)
  # } else if ( length(len) > 1 && sum(len < k) != length(len) ) {
  #   stop(paste("Some case regions have shorter size than k. I cannot decide",
  #              "on what to do."))
  # }

  # add buffer region
  case.region[, `:=`(start = start - control.relative.position[1],
                     end = end + control.relative.position[1])]

  # convert - to + strand
  case.region[strand == "-", strand := "+"]

  # trim out of range coordinates
  trimGenCoordinates("case.region", genome, remove = TRUE)

  # merge overlapped regions
  mergeGenCoordinate("case.region")

  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")

  # ----------------------------------------------------------------------------
  # SELECTION OF CONTROL REGIONS

  start.time <- Sys.time()

  cat("Selecting control regions...")

  # select control regions
  control.region <- env[[genomic.coordinate]][, .(
    start = c(start - control.relative.position[2], # [upstream] start
              end + control.relative.position[1]),  # (downstream) start
    end = c(start - control.relative.position[1],   # [upstream] end
            end + control.relative.position[2])     # (downstream) end
  ), by = .(chromosome, strand)]

  rm(genomic.coordinate, envir = env)
  gc()

  # convert "-" strand to "+" strand to merge them together
  control.region[strand == "-", strand := "+"]

  # trim out of range coordinates
  trimGenCoordinates("control.region", genome, remove = TRUE)

  # merge overlapped regions
  mergeGenCoordinate("control.region")

  # remove overlapping case-zone portion (taking too long - need to optimise)
  removeCaseZone("control.region", "case.region", strand.sensitive)

  # clear up the memory because we just memory copy genomic coordinates
  gc()

  # remove region with lower than size k
  control.region <- control.region[end - start + 1 >= k]

  # change to * strand to denote both strand
  control.region[, strand := "*"]

  # Save the coordinates
  fwrite(control.region, paste0(output, "/control_regions_coordinates.csv"))

  # remove unneeded object
  rm(case.region)
  gc()

  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")

  # ----------------------------------------------------------------------------
  # EXTRACTION OF CONTROL KMERS

  start.time <- Sys.time()

  cat("Extracting control kmers...")

  control.kmers <- extractKmers(control.region, genome, k, DNA.pattern)

  countReverseComplementKmers("control.kmers")

  if (!is.null(DNA.pattern)){

    # calculate expansion factor and pattern position
    expansion.factor <- (k-nchar(DNA.pattern[1]))/2
    pattern.pos <- seq(expansion.factor + 1,
                       expansion.factor + nchar(DNA.pattern[1]))

    control.kmers <- control.kmers[stri_sub(kmer, pattern.pos[1],
                                            pattern.pos[length(pattern.pos)
                                                        ]) %in% DNA.pattern]
  }

  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"), "\n")

  # ----------------------------------------------------------------------------
  # KMERS TABLE

  # rename column count to either case or control
  setnames(control.kmers, "count", "control")
  setnames(env[[case.kmers]], "count", "case")

  # merge the table
  kmers <- merge.data.table(control.kmers, env[[case.kmers]], all = TRUE)

  # change NA to 0
  kmers[is.na(control), control := 0]
  kmers[is.na(case), case := 0]

  # ----------------------------------------------------------------------------
  # MISC

  # restore back genomic coordinate table to original if changed
  if ("original_start" %in% colnames(env[[genomic.coordinate]])) {
    env[[genomic.coordinate]][, `:=`(start = original_start,
                                     end = original_end)]
    env[[genomic.coordinate]][, `:=`(original_start = NULL,
                                     original_end = NULL)]
  }

  return(kmers)
}
