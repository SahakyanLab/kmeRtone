## code to prepare `example_genome_coor` dataset goes here

library(data.table)
library(kmeRtone)

set.seed(1234)
example_genome_coor <- data.table(
    seqnames = "chr1",
    start = sample(
        x = 10000:100000, 
        size = 1000, 
        replace = FALSE
    ),
    width = 2
)

usethis::use_data(example_genome_coor, overwrite = TRUE)
