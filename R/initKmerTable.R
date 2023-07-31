#' Initialise k-mer table with all possible k-mers
#'
#' Initialise k-mer table with the following columns: kmer, pos_strand, and
#'    neg_strand
#'
#' @param k K-mer size. Limit to 15 because vector size is limited to
#'     .Machine$integer.max. For 9- to 15-mer, the kmer sequence is separated to
#'     two columns (kmer_part1 and kmer_part2) to reduce memory significantly.
#' @param central.pattern Central pattern(s) of the k-mer. Default is NULL.
#'
#' @return data.table with 3 columns: kmer, pos_strand, neg_strand
#'
#' @export
initKmerTable <- function (k, central.pattern=NULL, split.kmer=FALSE) {

  if (missing(k)) k <- NULL

  if (is.null(k)) {
    kmer.table <- data.table(kmer = character(),
                             pos_strand = numeric(),
                             neg_strand = numeric(), key = "kmer")
    if (split.kmer) {
      setnames(kmer.table, "kmer", "kmer_part1")
      kmer.table[, kmer_part2 := character()]
      setcolorder(kmer.table, c("kmer_part1", "kmer_part2"))
    }
  } else if (k > 15) {
    stop("k is limited to 15 because vector size is limited to ",
         ".Machine$integer.max")
  } else if (!split.kmer) {

    if (is.null(central.pattern)) {

      possible.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
      possible.kmers <- possible.kmers[, do.call(paste0, .SD)]
      kmer.table <- data.table(kmer = possible.kmers, pos_strand = 0,
                               neg_strand = 0, key = "kmer")

    } else {

      lens <- unique(nchar(central.pattern))

      kmer.table <- lapply(lens, function(len) {

        flank.size <- (k - len) / 2
        pattern.pos <- (flank.size + 1):(flank.size + len)

        flank.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")),
                                       flank.size))
        flank.kmers <- flank.kmers[, do.call(paste0, .SD)]
        flank.kmers <- do.call(CJ, rep(list(flank.kmers), 2))

        ptns <- central.pattern[nchar(central.pattern) == len]

        possible.kmers <-
          lapply(ptns, function(ptn) flank.kmers[, .(V1, ptn, V2)]) |>
          rbindlist()

        possible.kmers <- possible.kmers[, do.call(paste0,.SD)]
        kmer.table <- data.table(kmer = possible.kmers, pos_strand = 0,
                                 neg_strand = 0, key = "kmer")

        return(kmer.table)
      }) |> rbindlist()

      kmer.table <- unique(kmer.table)

    }

    # Separate k-mer sequence to two parts (sweet spot).
  } else if (split.kmer) {

    k1 <- k %/% 2
    k2 <- k - k1

    part1 <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k1))
    part1 <- part1[, do.call(paste0,.SD)]
    if (k2 == k1) {
      part2 <- part1
    } else {
      part2 <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k2))
      part2 <- part2[, do.call(paste0,.SD)]
    }

    kmer.table <- do.call(CJ, list(kmer_part1 = part1, kmer_part2 = part2))
    kmer.table[, c("pos_strand", "neg_strand") := 0]
    setkey(kmer.table, kmer_part1, kmer_part2)

    # Filter out different central pattern
    if(!is.null(central.pattern)) {

      lens <- unique(nchar(central.pattern))

      kmer.table <- lapply(lens, function(len) {

        flank.size <- (k - len) / 2
        pattern.pos <- (flank.size + 1):(flank.size + len)

        pattern.size.1 <- sum(pattern.pos <= k1)

        ptns <- central.pattern[nchar(central.pattern) == len]

        kmer.table <- lapply(ptns, function(ptn) {

          pattern.1 <- stri_sub(ptn, 1, pattern.size.1)
          pattern.2 <- stri_sub(ptn, pattern.size.1 + 1)

          kmer.table <- kmer.table[
            stri_detect_regex(kmer_part1, paste0(pattern.1, "$")) &
              stri_detect_regex(kmer_part2, paste0("^", pattern.2))]

          return(kmer.table)
        }) |> rbindlist()

        return(kmer.table)
      }) |> rbindlist()
    }
  }

  setkeyv(kmer.table, names(kmer.table)[grepl("kmer", names(kmer.table))])

  return(kmer.table)
}