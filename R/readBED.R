#' Read a BED file. One-based indexing is enforced.
#'
#' @param bed.path A path to a BED file.
#'
#' @return data.table.
#' 
#' @importFrom data.table setnames set fread
#'
#' @export
readBED <- function(bed.path) {

  # By convention, these are BedFile columns in a strict order.
  bed.cols <- c("chrom", "chromStart", "chromEnd", "name", "score",
                "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
                "blockSizes", "blockStarts")

  bed <- fread(bed.path, showProgress = FALSE)

  # Assign column names
  setnames(bed, colnames(bed), bed.cols[1:ncol(bed)])

  # Change to one-based indexing
  col.nums <- which(stri_detect_fixed(colnames(bed), "Start"))
  set(bed, j = col.nums, value = bed[, ..col.nums] + 1)

  return(bed)
}
