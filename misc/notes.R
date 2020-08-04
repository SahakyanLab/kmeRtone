

# initiate all possible kmers
# initiate count; should be rep(NA, all.possible.kmer)
kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
pattern.idx <- kmers[, do.call(paste0,.SD) %in% DNA.pattern,
                     .SDcols = pattern.pos]
kmers <- kmers[pattern.idx]
kmers <- kmers[, do.call(paste0,.SD)]
names(kmers) <- kmers
kmers[names(kmers)] <- 0
class(kmers) <- "numeric"



# To divide and process table by 100k rows later - mayble help with memory efficiency.
#    - Nothing changed for 13 million rows. Maybe above 13 million rows.
table.chunks <- gl(nrow(genomic.coordinate)/100000, 100000, nrow(genomic.coordinate))

# calculate expandion factor and pattern position
expansion.factor <- (k-nchar(DNA.pattern))/2
pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern))


# To-do list
# 1. Make extractKmers more memory efficient [DONE]
# 2. 