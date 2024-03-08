kmeRtone::kmeRtone(
    case.coor.path="../../data-raw/", 
    genome.name="hg19", 
    strand.sensitive=FALSE, 
    k=4,
    ctrl.rel.pos=c(80, 500),
    case.pattern=NULL,
    single.case.len=2,
    output.dir="./output",
    module="score",
    rm.case.kmer.overlaps=FALSE,
    merge.replicate=TRUE, 
    verbose=TRUE
)

ctrl_file <- data.table::fread("./output/score_4-mers.csv")
case_sum <- sum(ctrl_file$case)

testthat::test_that('case sum from scoring function', {
  testthat::skip_on_cran()
  testthat::expect_equal(193020, case_sum)
})

control_sum <- sum(ctrl_file$control)
testthat::test_that('control sum from scoring function', {
  testthat::skip_on_cran()
  testthat::expect_equal(3746116, control_sum)
})