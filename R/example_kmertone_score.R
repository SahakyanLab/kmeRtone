#' @title Example 2-mer enrichment/depletion scores
#' 
#' @description 
#' Below is an example code that generates random genomic coordinates 
#' and runs the default kmeRtone `SCORE` function to quantify the 
#' k-meric enrichment and depletion. 
#' 
#' @examples 
#' \dontrun{
#' library(data.table)
#' library(kmeRtone)
#' 
#' # 1. Randomly generate genomic positions and save results
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
#' # 2. Run kmeRtone `score` function
#' kmeRtone::kmeRtone(
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
#' }
#' 
#' @format A data frame with 1001 rows and 3 columns
#' \describe{
#'      \item{case}{Case k-mers, e.g. damage k-mer counts}
#'      \item{case_skew}{Case k-mers skews, e.g. skew of the damage k-mers counts}
#'      \item{control}{control k-mers, e.g. damage k-mer counts}
#'      \item{control_skew}{control k-mers skews, e.g. skew of the damage k-mers counts}
#'      \item{kmer}{K-meric sequence}
#'      \item{z}{Intrinsic susceptibility z-score for each k-mer}
#' }
#' @source \url{https://github.com/SahakyanLab/kmeRtone/blob/master/README.md}
#'
"example_kmeRtone_score"