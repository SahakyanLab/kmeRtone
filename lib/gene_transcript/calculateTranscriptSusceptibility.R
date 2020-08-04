calculateTranscriptSusceptibility <- function(element.coordinate, genome,
                                              ref.kmers, cutoff,
                                              output=NULL, verbose=FALSE,
                                              ncpu=1) {

  # AIM: Get kmers of each gene of element from sense and antisense strands.
  # This function relies on these columns: name, (optionally name2), element,
  # chromosome, strand, start, end.
  # It took 4 hours for NCBI RefSeq hg38 table (19,367 transcripts).

  # element.coordinate <data.table>  Genomic coordinate of transcript elements:
  #                                  CDS, intron, UTR, etc.
  # genome               <genome>    Genome-class object
  # ref.kmers          <data.table>  Kmers and respective score
  # score                <string>    Type of score to use/ score-column name
  # cutoff                <int>      Percentage of cutoff score to use.
  # output               <string>    [optional] Output filename.
  # verbose             <boolean>    Verbose message.

  # Dependencies
  #   Package   : data.table
  #   Functions : extractKmers, distributeChunk2

  if (cutoff > 100 | round(cutoff) < 0)
    stop("Cutoff should be between 0 and 100.")

  if (is.character(ref.kmers)) ref.kmers <- fread(ref.kmers)

  k <- ref.kmers[, unique(nchar(kmer))]

  if (element.coordinate[end - start + 1 >= k, .N] == 0)
    stop("No region longer than k.")

  # Sort the table
  columns <- colnames(element.coordinate)
  key.columns <- c("name", if("name2" %in% columns) "name2" else NULL,
                   "element")
  setkeyv(element.coordinate, c("chromosome", "strand", key.columns))

  # DEV - mark for deletion
  #assign("kmers.kmertone", NULL, envir = .GlobalEnv)

  # Main calculation helper function
  countSusceptibility <- function() {

    susceptibility <- element.coordinate[end - start + 1 >= k, {

      if (verbose) cat(name, "\t", element, "\n")

      # Extract kmers from the sense strand.
      kmers <- extractKmers(.SD, genome, k)

      # Total kmers
      total_kmers <- as.numeric(kmers[, sum(count)])

      # Include unresolved base N
      actual_total_kmers <- as.numeric(sum((end - start + 1) - k + 1))

      # Count antisense kmers
      countReverseComplementKmers("kmers", update.count = FALSE)

      setnames(kmers, c("count", "rc.count"), c("sense", "antisense"))

      # DEV - mark for deletion
      #kmers.kmertone[[element]] <<- kmers
      #print("Dah siapppp!")

      # Count susceptible kmers for every score
      scores <- colnames(ref.kmers)[4:ncol(ref.kmers)]

      for (score in scores) {

        # Total susceptible kmers in list format: list(sense, antisense)
        assign(paste0("susceptible_kmers_", score),
               as.numeric(kmers[kmer %in%
                            ref.kmers[!is.na(eval(parse(text = score)))]$kmer,
                          c(sum(sense), sum(antisense))]))

        setkeyv(ref.kmers, score)
        tail.percent <- (100 - cutoff) / 2
        high.risk.kmers <- ref.kmers[!is.na(eval(parse(text = score)))][
                          (round((100 - tail.percent) / 100 * .N) + 1):.N]$kmer
        low.risk.kmers <- ref.kmers[!is.na(eval(parse(text = score)))][
                              1:(round(tail.percent / 100 * .N) - 1)]$kmer

        # Total high-risk kmers in list format: c(sense, antisense)
        assign(paste0("high_risk_kmers_", score),
               as.numeric(kmers[kmer %in% high.risk.kmers,
                          c(sum(sense), sum(antisense))]))


        # Total low-risk kmers in vector format: c(sense, antisense)
        assign(paste0("low_risk_kmers_", score),
               as.numeric(kmers[kmer %in% low.risk.kmers,
                          c(sum(sense), sum(antisense))]))

      }

      column.names <- c("strand", "total_kmers", "actual_total_kmers",
                      CJ(c("susceptible_kmers_", "high_risk_kmers_",
                               "low_risk_kmers_"), scores)[, do.call(paste0, .SD)])
      strand <- c("sense", "antisense")

      sapply(column.names, function(n) eval(parse(text = n)))

    }, by = key.columns]

    return(susceptibility)
  }

  initial.thread <- getDTthreads()
  setDTthreads(1)

  time.start <- Sys.time()
  if (ncpu == 1) {

    susceptibility <- countSusceptibility()

  } else {

    chunk <- distributeChunk(element.coordinate[end - start + 1 >= k, .N], ncpu)
    dt.list <- lapply(seq_along(chunk$start), function(i) {
                 element.coordinate[end - start + 1 >= k
                                    ][chunk$start[i]:chunk$end[i]]
               })
    genome.list <- lapply(dt.list, function(dt) {
      genome <- genome[dt[, unique(chromosome)]]
      class(genome) <- "genome"
      return(genome)
    })

    to.export.variable <- c("key.columns", "ref.kmers", "cutoff",
                            "k", "verbose")
    to.export.function <- c("countSusceptibility", "extractKmers",
                            "distributeChunk2", "countReverseComplementKmers",
                            "reverseComplement")

    susceptibility <-
      foreach(element.coordinate = dt.list, genome = genome.list,
              .combine = "rbind", .packages = c("data.table", "stringi"),
              .export = c(to.export.variable, to.export.function),
              .noexport = ls()[!ls() %in% to.export.variable]) %dopar% {

      setDTthreads(1)
      susceptibility <- countSusceptibility()

      return(susceptibility)
    }
    # Some of the multiple element block could be separated.
    if (sum(duplicated(susceptibility[, 1:3])) > 0) {

      susceptibility <-
        susceptibility[, lapply(.SD, sum), by = c(key.columns, "strand")]
    }
  }

  time.diff <- Sys.time() - time.start
  cat("Time taken\t:", time.diff, attr(time.diff, "units"), "\n")

  if (!is.null(output))
    fwrite(susceptibility, paste0(output, "/", "transcript_susceptibility.csv"))

  setDTthreads(initial.thread)

  return(susceptibility)
}
