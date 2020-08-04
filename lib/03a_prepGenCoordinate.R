prepGenCoordinate <- function(genomic.coordinate, strand.sensitive, genome,
                              env=parent.frame()) {

  # Dependencies:
  #     Kmertone variables : genomic.coordinate, strand.sensitive
  #     Package            : data.table

  essential.columns <- c("chromosome", "start", "end", "strand")

  # processing a string format
  if (class(env[[genomic.coordinate]])[1] == "character") {

    filename <- env[[genomic.coordinate]]

    # genomic.coordinate can be replicates
    env[[genomic.coordinate]] <- lapply(env[[genomic.coordinate]], function(i) {
                                        fread(i, showProgress = FALSE)})

    # BED files format
    if (sum(grepl("\\.bed", filename)) == length(filename)) {

      cat("--BED file format detected from the filename extension.\n")

      for (i in seq_along(env[[genomic.coordinate]])) {

        if (ncol(env[[genomic.coordinate]][[i]]) >= 4) {
          # assume strand sensitive i.e. has strand information
          # [TAGGED] Can do a better detection i.e. directly check for "+" "-"

          if (!strand.sensitive) {
            warning(paste0("--The strand mode specified is sensitive. Strand",
                           "information will be ignored."))
          }

          for (colnum in 4:(ncol(env[[genomic.coordinate]][[i]]))) {

            element <- unique(env[[genomic.coordinate]][[i]][[colnum]])
            if (sum((element %in% c("+", "-", "*"))) == length(element)) {
              strand.column <- paste0("V", colnum)}

          }
          # Only retain four essential columns.
          if (ncol(env[[genomic.coordinate]][[i]]) > 4) {
            env[[genomic.coordinate]][[i]][,
               colnames(env[[genomic.coordinate]][[i]])[
                  !colnames(env[[genomic.coordinate]][[i]]) %in%
                       c("V1", "V2", "V3", strand.column)] := NULL]
          }

          setnames(env[[genomic.coordinate]][[i]], c("V1", "V2", "V3",
                                                     strand.column),
                   essential.columns)

          # bed use zero-based index and open-end index e.g. [0,10) which means
          # from index 0 until index 9
          # change to R indexing
          env[[genomic.coordinate]][[i]][, start := start + 1]

        } else {

          if (strand.sensitive) {
            stop(paste0("You specify strand sensitive but there is no strand",
                 "information in the table."))
          } else if (!strand.sensitive) {

            # strand insensitive i.e. has no strand information
            env[[genomic.coordinate]][[i]][, colnames(env[[
                    genomic.coordinate]][[i]])[!colnames(env[[
                       genomic.coordinate]][[i]]) %in%
                            c("V1", "V2", "V3")] := NULL]

            setnames(env[[genomic.coordinate]][[i]], c("V1", "V2", "V3"),
                     essential.columns[-4])
            env[[genomic.coordinate]][[i]][, start := start + 1]
            env[[genomic.coordinate]][[i]][, strand := "*"]
          }
        }
      }
    }
  }

  if (class(env[[genomic.coordinate]])[1] == "list" &
      length(env[[genomic.coordinate]]) == 1){
    env[[genomic.coordinate]] <- env[[genomic.coordinate]][[1]]
    gc()
  }

  # if R object list is given/generated, assume a list of data.table
  if (class(env[[genomic.coordinate]])[1] == "list") {

    cat("--Detecting", length(env[[genomic.coordinate]]), "replicates\n")

    # rename columns
    for (i in seq_along(env[[genomic.coordinate]])) {

      columns.in.table <- colnames(env[[genomic.coordinate]][[i]])

      # If essential columns not exist
      if(sum(essential.columns %in% columns.in.table) < 3){

        if (ncol(env[[genomic.coordinate]][[i]]) > 3) {
          setnames(env[[genomic.coordinate]][[i]], columns.in.table[1:4],
                   essential.columns)
        } else {
          setnames(env[[genomic.coordinate]][[i]], columns.in.table[1:3],
                   essential.columns[-4])
          env[[genomic.coordinate]][[i]][, strand := "*"]
        }
      }

      non.essential.columns <- columns.in.table[!columns.in.table %in%
                                                essential.columns]

      # Remove other than essential columns
      if(!sum(essential.columns %in% columns.in.table) %in% 3:4){
        env[[genomic.coordinate]][[i]][, eval(non.essential.columns) := NULL]
      }

    }

    # add column replicate_number if not exist yet
    if(sum(grepl("replicate_", colnames(env[[genomic.coordinate]]))) < 1) {

      cnt <- 1
      for (i in 1:length(env[[genomic.coordinate]])) {
        env[[genomic.coordinate]][[i]][, paste0("replicate_", cnt) := 1]
        cnt <- cnt + 1
      }

    }

    essential.columns <- essential.columns[columns.in.table %in%
                                           essential.columns]

    cat("--Merging the genomic coordinate tables...\n")
    dt <- Reduce(function(dt1, dt2) merge.data.table(dt1[, ], dt2,
                                           by = essential.columns, all = TRUE),
                 x = env[[genomic.coordinate]])

    #fwrite(dt, "data/consolidated_table.csv")

    # reassign the consolidated table to genomic.coordinate
    env[[genomic.coordinate]] <- dt

    # memory copy, so gc() to clear memory
    gc()

  } else {

    if (ncol(env[[genomic.coordinate]]) > 3) {
      setnames(x = env[[genomic.coordinate]],
               old = colnames(env[[genomic.coordinate]])[1:4],
               new = c("chromosome", "start", "end", "strand"))
    } else if (ncol(env[[genomic.coordinate]]) == 3) {
      setnames(x = env[[genomic.coordinate]],
               old = colnames(env[[genomic.coordinate]])[1:3],
               new = c("chromosome", "start", "end"))
      env[[genomic.coordinate]][, strand := "*"]
    }


  }

  setkey(env[[genomic.coordinate]], chromosome)
  if (sum(!env[[genomic.coordinate]][, unique(chromosome)] %in%
          names(genome)) > 0) {

    message(paste0("Genome chromosome names: ", paste(names(genome),
                                                      collapse = ", ")))
    message(paste0("Genomic coordinate table chromosome names: ",
                   paste(env[[genomic.coordinate]][, unique(chromosome)],
                         collapse = " ")))
    stop(paste0("Chromosome names are not consistent between genome and",
                "genomic coordinate table"))

  }
}
