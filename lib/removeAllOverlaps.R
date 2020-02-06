removeAllOverlaps <- function(genomic.coordinate, env) {
  # Remove all overlapping regions. This is used for removing
  # overlapping case regions.
  
  # genomic.coordinate    <string>        A variable name pointing to data.table of
  #                                       genomic coordinate.
  # env                  <environment>    An environment object where the table exists.
  
  ## locate continuous overlapping regions --> group them together and assign group number
  setkey(env[[genomic.coordinate]], chromosome, strand, start, end)
  env[[genomic.coordinate]][, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < 0)),
                         by = list(chromosome, strand)]
  
  ## save overlaps
  fwrite(env[[genomic.coordinate]][, if (.N > 1) .(start, end),
                                   by = .(chromosome, strand, group)],
         "data/overlapping_case_region.csv")
  
  ## remove overlaps
  env[[genomic.coordinate]] <- env[[genomic.coordinate]][, if (.N == 1) .(start, end),
                                                   by = .(chromosome, strand, group)]
  # memory copy, so gc()
  gc()
  
  # remove the columns
  env[[genomic.coordinate]][, group := NULL]
  setcolorder(env[[genomic.coordinate]], c("chromosome", "start", "end", "strand"))
  
}