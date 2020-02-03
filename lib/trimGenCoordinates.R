trimGenCoordinates <- function(genomic.coordinate, genome, env=parent.env(), omit.NA=FALSE) {
  # DESCRIPTION
  # Trim out of range genomic coordinates
  #
  # genomic.coordinate    <string>     A variable of data.table; must have these columns: chromosome, start, end, strand
  # genome                <string>     A variable of genome
  # env                 <environment>  An environment where genomic.coordinate and genome located
  
  # remove -ve ranges on both start and end
  env[[genomic.coordinate]][start <= 0 & end <= 0, c("start", "end") := NA]

  # trim out of range genomic coordinate
  # less than genome chromosome start number i.e. < 1
  env[[genomic.coordinate]][start < 1, start := 1]
  env[[genomic.coordinate]][end < 1, end := 1]
  
  # more than chromosome end number
  env[[genomic.coordinate]][start > nchar(env[[genome]][[chromosome]]) & end > nchar(env[[genome]][[chromosome]]),
                            c("start", "end") := NA]
  env[[genomic.coordinate]][start > length(env[[genome]][[chromosome]]), start := length(env[[genome]][[chromosome]])]
  env[[genomic.coordinate]][end > length(env[[genome]][[chromosome]]), end := length(env[[genome]][[chromosome]])]
  
  if (omit.NA == TRUE) {
    genomic.coordinate <- na.omit(env[[genomic.coordinate]])
    return(genomic.coordinate)
  }
}
