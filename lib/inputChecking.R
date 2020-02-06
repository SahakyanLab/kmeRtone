inputChecking <- function(env=kmertone.env) {
  
  # genomic.coordinate is checked in prepGenCoordinate function.
  # genome.name, genome.path, genome.prefix, genome.suffix are checked in prepGenome function.
  
  # Dependencies
  #    Kmertone variables: DNA.pattern, k, control.relative.position, strand.mode, ncpu
  
  if ((!class(env$DNA.pattern) %in% c("character", "NULL")) | (!length(env$DNA.pattern) %in% 0:1)) {
    stop("Please input a single DNA pattern in a string format.")
  }
  
  if (!class(env$k) %in% c("numeric", "integer") || env$k%%1 != 0) {
    stop("Please input a round number for k.")
  }
  
  if (length(env$control.relative.position) != 2 ||
      sum(!class(env$control.relative.position) %in% c("numeric", "integer")) > 0 ||
      sum(env$control.relative.position %% 1 != 0) > 0) {
    stop("Please input a control relative position in a vector of 2 round numbers.")
  }
  
  if (sum(class(env$strand.mode) != "character") > 0 ||
      sum(!env$strand.mode %in% c("sensitive", "insensitive")) > 0) {
    stop("Please strand mode in either sensitive or insensitive")
  }
  
  if (!class(env$ncpu) %in% c("numeric", "integer") || length(env$ncpu) > 1 || env$ncpu%%1 != 0) {
    stop("Please input correct number of cpu to use.")
  }
}