#' Write a BED file. Zero-based indexing is enforced.
#'
#' @param bed A BED `data.table`.
#' @param output.filename An output BED filename.
#'
#' @export
readBED <- function(bed, output.filename) {

  # By convention, these are BedFile columns in a strict order.
  bed.cols <- c("chrom", "chromStart", "chromEnd", "name", "score",
                "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
                "blockSizes", "blockStarts")

  # Check required columns: chrom, chromStart, and chromEnd
  if (all(bed.cols[1:3] %in% colnames(bed))) {
    setcolorder(bed, bed.cols[1:3])
  } else {
    if (!(is.character(bed[1][[1]]) & is.numeric(bed[1][[2]]) &
      is.numeric(bed[1][[3]]))) {
        stop(paste0("Failed to detect required columns: chrom, chromStart, ",
                    "and chromEnd."))
      }
  }

   # Change back to zero-based indexing
  col.nums <- which(stri_detect_regex(colnames(bed), "[Ss]tart"))
  set(bed, j = col.nums, value = bed[, ..col.nums] + 1)

  bed <- fwrite(output.filename, col.names = FALSE, showProgress = FALSE)

  # Change back to zero-based indexing
  set(bed, j = col.nums, value = bed[, ..col.nums] - 1)

}
