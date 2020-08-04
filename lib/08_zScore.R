zScore <- function(kmers, DNA.pattern, env=parent.frame()) {
  # Calculate Z score and update the table
  # Claudis's thesis equation 1
  # Different pattern will be calculated separately

  # kmers   <string>   A variable name pointing to kmers<data.table>. The

  # Dependency
  #    Packages: data.table, stringi

  env[[kmers]][, z := {

    # total case count (n)
    total.case <- sum(case)

    # proportion control (p)
    p.control <- control / sum(control)

    # predicted case distribution (np)
    case.predict <- total.case * p.control

    z <- (case - case.predict) / sqrt( case.predict * (1 - p.control) )

    z
  }]

  # Length of kmer
  k <- env[[kmers]][, unique(nchar(kmer))]
  expansion.factor <- (k - unique(nchar(DNA.pattern))) / 2

  if (length(DNA.pattern) > 1) {

    # Add column DNA_pattern
    env[[kmers]][, DNA_pattern := stri_sub(kmer, expansion.factor + 1,
                                           k - expansion.factor)]

    for (dna.pattern in DNA.pattern) {

      env[[kmers]][DNA_pattern == dna.pattern, paste0("z_", dna.pattern) := {

        # total case count (n)
        total.case <- sum(case)

        # proportion control (p)
        p.control <- control / sum(control)

        # predicted case distribution (np)
        case.predict <- total.case * p.control

        z <- (case - case.predict) / sqrt( case.predict * (1 - p.control))

        z
      }]
    }

    env[[kmers]][, DNA_pattern := NULL]

  }
}
