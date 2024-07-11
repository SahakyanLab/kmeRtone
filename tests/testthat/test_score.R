library(data.table)
library(kmeRtone)
temp_dir <- tempdir()

set.seed(1234)
temp_files <- character(22)
for(chr in 1:22){
    genomic_coor <- data.table(
        seqnames = paste0("chr", chr),
        start = sample(
            x = 10000:100000, 
            size = 1000, 
            replace = FALSE
        ),
        width = 2
    )

    f <- file.path(temp_dir, paste0("chr", chr, ".csv"))
    fwrite(genomic_coor, f)
    temp_files[chr] <- f
}

# run kmertone score function
kmeRtone::kmeRtone(
    case.coor.path = temp_dir, 
    genome.name = "hg19", 
    strand.sensitive = FALSE, 
    k = 2,
    ctrl.rel.pos = c(80, 500),
    case.pattern = NULL,
    single.case.len = 2,
    output.dir = temp_dir,
    module = "score",
    rm.case.kmer.overlaps = FALSE,
    merge.replicate = TRUE, 
    kmer.table = NULL,
    verbose = TRUE
)

ctrl_file <- data.table::fread(paste0(temp_dir, "/score_2-mers.csv"))
case_sum <- sum(ctrl_file$case)

testthat::test_that('case sum from scoring function', {
  testthat::skip_on_cran()
  testthat::expect_equal(24970, case_sum)
})

control_sum <- sum(ctrl_file$control)
testthat::test_that('control sum from scoring function', {
  testthat::skip_on_cran()
  testthat::expect_equal(385028, control_sum)
})

rm_files <- file.remove(temp_files)