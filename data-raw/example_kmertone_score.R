## code to prepare `example_kmeRtone_score` dataset goes here

example_kmeRtone_score <- read.csv("data-raw/score_2-mers.csv")

usethis::use_data(example_kmeRtone_score, overwrite = TRUE)