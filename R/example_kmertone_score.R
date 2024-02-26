#' @title Example 2-mer enrichment/depletion scores
#' 
#' @description 
#' Below is an example code that generates random genomic coordinates 
#' and runs the default kmertone `SCORE` function to quantify the 
#' k-meric enrichment and depletion. 
#' 
#' library(data.table)
#' library(kmertone)
#' 
#' 1. Randomly generate genomic positions and save results
#' dir.create("./data", showWarnings = FALSE)
#' 
#' set.seed(1234)
#' for(chr in 1:22){
#'     genomic_coor <- data.table::data.table(
#'         seqnames = paste0("chr", chr),
#'         start = sample(
#'             x = 10000:10000000, 
#'             size = 100000, 
#'             replace = FALSE
#'         ),
#'         width = 2
#'     )
#' 
#'     data.table::fwrite(
#'         genomic_coor, 
#'         paste0("./data/chr", chr, ".csv")
#'     )
#' }
#' 
#' #' 2. Run kmertone `score` function
#' kmertone::kmertone(
#'     case.coor.path="./data", 
#'     genome.name="hg19", 
#'     strand.sensitive=FALSE, 
#'     k=2,
#'     ctrl.rel.pos=c(80, 500),
#'     case.pattern=NULL,
#'     single.case.len=2,
#'     output.dir="output",
#'     module="score",
#'     rm.case.kmer.overlaps=FALSE,
#'     merge.replicate=TRUE, 
#'     kmer.table=NULL,
#'     verbose=TRUE
#' )
#' 
#' @format A data frame with 1001 rows and 3 columns
#' \describe{
#'   \item{seqnames}{Chromosome number of the recorded biological event, e.g. DNA strand breaks}
#'   \item{start}{5' start position of the recorded biological event}
#'   \item{width}{Sequence width of the recorded biological event, e.g. 2 for a DNA strand break}
#' }
#' @source \url{https://github.com/SahakyanLab/kmertone/tree/master/README.md}
#'
"example_kmertone_score"