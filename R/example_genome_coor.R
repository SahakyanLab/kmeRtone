#' @title Example genome coordinate file
#' 
#' @description 
#' Below is an example code that generates random genomic coordinates.
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
#' for(chr in 1){
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
#' }
#' 
#' @format A data frame with 1001 rows and 3 columns
#' \describe{
#'   \item{seqnames}{Chromosome number of the recorded biological event, e.g. DNA strand breaks}
#'   \item{start}{5' start position of the recorded biological event}
#'   \item{width}{Sequence width of the recorded biological event, e.g. 2 for a DNA strand break}
#' }
#'
"example_genome_coor"