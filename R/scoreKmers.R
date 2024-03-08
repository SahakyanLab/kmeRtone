#' Function calculates the Z-score for each k-mer based on the observed 
#'  case counts and expected case counts under the null hypothesis.
#' 
#' @param kmer.table A data.table containing k-mer counts, where each row represents a 
#'  k-mer and columns "case" and "control" represent the counts in case and control samples respectively.
#' 
#' @return A modified version of the input `kmer.table` with an additional column 
#'  "z" containing the calculated Z-scores for each k-mer.
#' 
#' @export
scoreKmers <- function(kmer.table) {

  kmer.table[, z := {

    # total case count (n)
    total.case <- sum(case)

    # proportion control (p)
    p.control <- control / sum(control)

    # predicted case distribution (np)
    case.predict <- total.case * p.control

    z <- (case - case.predict) / sqrt( case.predict * (1 - p.control) )

    z
  }]

  return(kmer.table)
}