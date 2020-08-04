inputChecking <- function(DNA.pattern, k, control.relative.position,
                          strand.sensitive) {

  # genomic.coordinate is checked in prepGenCoordinate function.
  # genome.name, genome.path, genome.prefix, genome.suffix are checked in
  # prepGenome function.

  # Dependencies
  #    Kmertone variables: DNA.pattern, k, control.relative.position,
  #                        strand.sensitive, ncpu

  if ((!class(DNA.pattern) %in% c("character", "NULL"))) {
    stop("Please input DNA pattern in a string format.")
  }

  if (!class(k) %in% c("numeric", "integer") || k%%1 != 0) {
    if(k != "optimum") stop("Please input a round number for k.")
  }

  if (length(control.relative.position) != 2 ||
      sum(!class(control.relative.position) %in% c("numeric", "integer")) > 0 ||
      sum(control.relative.position %% 1 != 0) > 0) {
    stop("Please input a control relative position in a vector of 2 round
         numbers.")
  }

  if (sum(class(strand.sensitive) != "logical") > 0) {
      stop("Please input TRUE of FALSE for strand.sensitive")
  }

  if (!class(ncpu) %in% c("numeric", "integer") || length(ncpu) > 1 || 
      ncpu%%1 != 0) {
    stop("Please input correct number of cpu to use.")
  }
}
