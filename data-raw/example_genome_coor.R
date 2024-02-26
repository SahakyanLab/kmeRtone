## code to prepare `example_genome_coor` dataset goes here

example_genome_coor <- read.csv("data-raw/chr1.csv")

usethis::use_data(example_genome_coor, overwrite = TRUE)