#' @title Example 2-mer enrichment/depletion scores
#' 
#' @description 
#' Below is an example code that generates random genomic coordinates 
#' and runs the default kmeRtone `SCORE` function to quantify the 
#' k-meric enrichment and depletion. 
#' 
#' @examples 
#' \donttest{
#' # 1. Randomly generate genomic positions and save results
#' library(data.table)
#' library(kmeRtone)
#' temp_dir <- tempdir()
#' 
#' set.seed(1234)
#' temp_files <- character(1)
#' for(chr in 1){
#'     genomic_coor <- data.table(
#'         seqnames = paste0("chr", chr),
#'         start = sample(
#'             x = 10000:10000000, 
#'             size = 100000, 
#'             replace = FALSE
#'         ),
#'         width = 2
#'     )
#' 
#'     f <- file.path(temp_dir, paste0("chr", chr, ".csv"))
#'     fwrite(genomic_coor, f)
#'     temp_files[chr] <- f
#' }
#' 
#' # 2. Run kmeRtone score function
#' temp_dir_genome <- tempdir()
#' kmeRtone::kmeRtone(
#'     case.coor.path = temp_dir, 
#'     genome.name = "hg19", 
#'     genome.path = temp_dir_genome,
#'     strand.sensitive = FALSE, 
#'     k = 2,
#'     ctrl.rel.pos = c(80, 500),
#'     case.pattern = NULL,
#'     single.case.len = 2,
#'     output.dir = temp_dir,
#'     module = "score",
#'     rm.case.kmer.overlaps = FALSE,
#'     merge.replicate = TRUE, 
#'     kmer.table = NULL,
#'     verbose = TRUE
#' )
#' 
#' # 3. Clean up temporary files
#' rm_files <- file.remove(temp_files)
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