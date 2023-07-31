#' Build Coordinate object.
#'
#' The Coordinate object is capable of loading genomic coordinates on demand.
#'     Chromosome-specific coordinates can be called in a bracket
#'     e.g. coor["chr1"] or coor[c("chr1", "chr2")]. The coordinates can also
#'     be expanded to k-mer size equally on both flanks
#'     e.g. coor["chr1", state = "kmer", k = 9]
#'
#'
#' @param root.path A path to a directory containing either:
#'     (1) chromosome-separated coordinate files
#'         (multiple replicates is assumed for sub-folder) or
#'     (2) bedfile (multiple replicates is assumed for separate bedfiles).
#' @param single.len Single case length relevant when all coordinates have the
#'     same length. This is for memory optimization. Default is NULL.
#' @param is.strand.sensitive A boolean whether strand polarity matters.
#'     Default is TRUE.
#' @param merge.replicates Merge coordinate from different replicates.
#'     Default is TRUE. If not merging, duplicates will give weight to the
#'     k-mer counting. If add.col.rep, merged coordinate will contain
#'     column replicate e.g. "rep1&rep2".
#' @param rm.dup Remove duplicates in each replicate. Default is TRUE.
#' @param add.col.rep Add column replicate to the coordinate table.
#' @param is.kmer Is the coordinate refers to k-mer i.e. expanded case?
#'     Default is FALSE.
#' @param k Length of k-mer relevant only when is.kmer is TRUE.
#' @param load.limit Maximum number of coordinate data.table loaded on RAM.
#'     Default is 1.
#' @return Coordinate object.
#'
#' @examples
#' coor <- loadCoordinate(...)
#' coor["chr1"]
#' coor[c("chr1", "chr2")]
#' coor["chr1", state = "kmer", k = 8]
#'
#' @export
loadCoordinate <- function(root.path = NULL,
                           single.len = NULL,
                           is.strand.sensitive = TRUE,
                           merge.replicates = TRUE,
                           rm.dup = TRUE,
                           add.col.rep = FALSE,
                           is.kmer = FALSE,
                           k = NA,
                           ori.first.index = 1,
                           load.limit = 1) {

  if (!is.kmer) {
    k <- NA
  } else if (is.kmer & !is.numeric(k)) {
    stop("Please input numeric k.")
  }

  coor <- Coordinate$new(root.path, single.len, is.strand.sensitive,
                         merge.replicates, rm.dup, add.col.rep, is.kmer, k,
                         ori.first.index, load.limit)



  return(coor)
}
