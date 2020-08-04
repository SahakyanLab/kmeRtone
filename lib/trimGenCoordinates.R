trimGenCoordinates <- function(genomic.coordinate, genome, remove=FALSE,
                               env=parent.frame()) {
  # DESCRIPTION
  # Trim out of range genomic coordinates
  #
  # genomic.coordinate    <string>     A variable of data.table; must have these columns:
  #                                    chromosome, start, end, strand.
  # genome                <string>     A <genome> object.
  # env                 <environment>  An environment where genomic.coordinate and genome located
  # remove                <string>     If TRUE, out of range start AND end coordinate will be removed
  #                                    from the table. If FALSE, NA will be introduced in place of
  #                                    coordinate.
  
  if (class(genomic.coordinate)[1] != "character" | length(genomic.coordinate) != 1) {
    stop("Please input genomic.coordinate and genome as a variable name in string format.")
  }

  # remove -ve ranges on both start and end
  env[[genomic.coordinate]][start <= 0 & end <= 0, c("start", "end") := NA]
  
  # trim out of range genomic coordinate
  # less than genome chromosome start number i.e. < 1
  env[[genomic.coordinate]][start < 1, start := 1]
  env[[genomic.coordinate]][end < 1, end := 1]
  # more than chromosome end number
  env[[genomic.coordinate]][start > attr(genome, "length")[chromosome] & 
                              end > attr(genome, "length")[chromosome],
                            c("start", "end") := NA]
  env[[genomic.coordinate]][start > attr(genome, "length")[chromosome],
                            start := attr(genome, "length")[chromosome]]
  env[[genomic.coordinate]][end > attr(genome, "length")[chromosome],
                            end := attr(genome, "length")[chromosome]]
  
  if (remove == TRUE) {
    env[[genomic.coordinate]] <- na.omit(env[[genomic.coordinate]])
  }
}
