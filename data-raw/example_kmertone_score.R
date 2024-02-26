## code to prepare `example_kmertone_score` dataset goes here

example_kmertone_score <- read.csv("data-raw/score_2-mers.csv")

usethis::use_data(example_kmertone_score, overwrite = TRUE)