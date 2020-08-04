calculateOptimumK <- function(genomic.coordinate, genome, set.k, DNA.pattern,
                              strand.sensitive, env=parent.frame()) {
  # To calculate optimum k value. The calculation is based on
  # Claudia's thesis equation 3 and 4.
  # This is a stand-alone function - separate from kmertone function.

  # Dependencies:
  #     Packages   : data.table, stringi
  #     Functions  : prepGenome, prepGenCoordinate, addColumnSequence,
  #                  filterTable, getCaseKmers
  #     Information: proportion of precalculated A, C, G, and T of the whole
  #                  genome (hg19)

  DNA.pattern.size <- unique(nchar(DNA.pattern))
  set.k <- seq(DNA.pattern.size, 14, 2)

  kmers.list <- lapply(set.k, function(k) {

    cat("\nk = ", k, "...")
    kmers <- getCaseKmers("genomic.coordinate", genome, k, DNA.pattern,
                          strand.sensitive, remove.overlaps = TRUE, env)
    env[["genomic.coordinate"]][, c("start", "end") := list(original_start,
                                                            original_end)]

    return(kmers)
  })

  names(kmers.list) <- set.k

  # for (k in set.k) {
  #   fwrite(kmers.list[[as.character(k)]], paste0("data/case-kmers_k-", k,
  #                                                ".csv"))
  # }

  cat("[5] Calculating optimal k...\n")

  q <- sapply(set.k[-1], function(k) {

    q <- kmers.list[[as.character(k)]][, {

      total.case <- sum(count)

      nt.init <- stri_sub(kmer, 1, 1)
      seq.mid <- stri_sub(kmer, 2, k - 1)
      nt.last <- stri_sub(kmer, k, k)

      # take proportion from previous k
      p.mid <- kmers.list[[as.character(k - 2)]][, {

        names(count) <- kmer
        count[seq.mid] / sum(count)

      }]

      # proportion of precalculated A, C, G, and T of the whole genome (hg19)
      p.nt <- c(A = 0.295, T = 0.295, C = 0.205, G = 0.205)

      p.init <- p.nt[nt.init]
      p.last <- p.nt[nt.last]

      p = p.init * p.mid * p.last

      # equation 3 - to compute chi square
      Q = total.case * sum((count / total.case - p)**2 / p)

      q = (Q - (4**(k - 2) - 1)) / sqrt(2 * (4**(k - 2) - 1))

      q
    }]

    names(q) <- as.character(k)

    return(q)
  })

}
