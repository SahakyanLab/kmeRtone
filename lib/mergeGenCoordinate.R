mergeGenCoordinate <- function(genomic.coordinate, env=parent.frame()) {
  # must have chromosome, start, end, and strand columns
  # exact output as reduce function from GenomicRanges

  if (class(genomic.coordinate)[1] != "character" | length(genomic.coordinate) != 1) {
    stop("Please input genomic.coordinate as a variable name in string format.")
  }
  
  # sort the table
  setkey(env[[genomic.coordinate]], chromosome, strand, start, end)

  # locate the overlapping
  env[[genomic.coordinate]][, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < -1)),
                             by = list(chromosome, strand)]

  # merge the overlaps
  env[[genomic.coordinate]] <- env[[genomic.coordinate]][, .(start = min(start), end = max(end)),
                                                           by = list(chromosome, strand, group)]

  # remove column group
  env[[genomic.coordinate]][, group := NULL]
  
  # rearrange the columns
  setcolorder(env[[genomic.coordinate]], c("chromosome", "start", "end", "strand"))

}