

# initiate all possible kmers
# initiate count; should be rep(NA, all.possible.kmer)
case.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), 10))
pattern.idx <- case.kmers[, rowSums(.SD == strsplit(DNA.pattern, "")[[1]]) == nchar(DNA.pattern),
                          .SDcols = pattern.pos]
case.kmers <- case.kmers[pattern.idx]
case.kmers <- case.kmers[, paste(.SD, collapse = ""), by = 1:nrow(case.kmers)][[2]]
names(case.kmers) <- case.kmers
case.kmers[names(case.kmers)] <- 0
class(case.kmers) <- "numeric"



# To divide and process table by 100k rows later - mayble help with memory efficiency.
#    - Nothing changed for 13 million rows. Maybe above 13 million rows.
table.chunks <- gl(nrow(genomic.coordinate)/100000, 100000, nrow(genomic.coordinate))

# calculate expandion factor and pattern position
expansion.factor <- (k-nchar(DNA.pattern))/2
pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern))