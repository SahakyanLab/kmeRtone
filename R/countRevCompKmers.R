#' Count reverse complement sequence from its opposite strand.
#' Build for k-mer table generated from initKmerTable function but applicable to
#' others with the same format.
#'
#' @param kmer.table A data.table of k-mer with at least 3 columns: kmer,
#'     pos_strand, and neg_strand. Splitted k-mer columns: kmer_part1 and
#'     kmer_part2 is supported.
#'
#' @return Updated k-mer table.
#'
#' @importFrom data.table setkey
#' @importFrom stringi stri_sub
#' @importFrom Biostrings reverseComplement
#' 
#' @export
countRevCompKmers <- function(kmer.table) {

  # Error checking. Count in neg_strand should be rationally zero
  if (kmer.table[, sum(neg_strand)] > 0) {
    stop("Count on negative strand is already populated!")
  } else {
    kmer.table[, neg_strand := NULL]
  }

  if (any(grepl("kmer_part", names(kmer.table)))) {

    setkey(kmer.table, kmer_part1, kmer_part2)

    # Only suitable for k < 14
    # Find reverse complement
    rc <- kmer.table[, .(reverseComplement(kmer_part2),
                         reverseComplement(kmer_part1),
                         pos_strand)]
    setnames(rc, names(rc), c("kmer_part1", "kmer_part2", "neg_strand"))

    # By kmer.table convention, the first part of kmer is shorter if k is odd
    if (rc[1, nchar(kmer_part1) != nchar(kmer_part2)]) {
      rc[, `:=`(kmer_part2 = paste0(stri_sub(kmer_part1,
                                             k1 <- nchar(kmer_part1[1]), k1),
                                    kmer_part2),
                kmer_part1 = stri_sub(kmer_part1, 1, k1 - 1)
      )]
    }
    setkey(rc, kmer_part1, kmer_part2)

    # Add new kmer from rc.seq. Memory copy unfortunately...
    if (kmer.table[!rc, .N] > 0) {
      kmer.table <- merge(kmer.table, rc, all = T)
      setkey(kmer.table, kmer_part1, kmer_part2)

      # Change NA to 0 for non-existing kmer
      kmer.table[is.na(pos_strand), pos_strand := 0]
      kmer.table[is.na(neg_strand), neg_strand := 0]
    } else {
      kmer.table[, neg_strand := rc$neg_strand]
    }

    return(kmer.table)
  }

  # Find reverse complement
  kmer.table[, rc.seq := reverseComplement(kmer)]

  # Some reverse complement sequence may not exist in the plus strand, so add
  # the new sequence if any
  if(kmer.table[!rc.seq %in% kmer, .N > 0]){
    kmer.table <- rbind(kmer.table,
                        kmer.table[!rc.seq %in% kmer,
                                   .(kmer = rc.seq, pos_strand = 0,
                                     rc.seq = kmer)])
  }

  # Count complementary sequence
  rc.count <- kmer.table[match(rc.seq, kmer), pos_strand]

  # Because some reverse complement sequence do not exist in the plus strand,
  # match() results in NAs. Change NA to zero
  if (any(is.na(rc.count))) {
    rc.count[is.na(rc.count)] <- 0
  }

  # Add column rc.count
  kmer.table[, neg_strand := rc.count]

  # Remove rc.seq column
  kmer.table[, rc.seq := NULL]

  return(kmer.table)
}
