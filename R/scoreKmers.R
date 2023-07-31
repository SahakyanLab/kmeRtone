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