#' Convert a BED file to chromosome-separated csv files.
#'
#' @param bed.path A path to a BED file.
#' @param output.path Output directory path. It should be empty directory.
#'     Default to coordinate/
#'
#' @return None
#'
#' @export
bedToCoor <- function(bed.path, output.path="coordinate/", compress = TRUE) {

  # Convert bedfile to normal 1-index csv file, separated by chromosome in a
  # a folder.

  if (length(list.files(output.path)) > 0) stop(output.path, " is not empty.")
  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

  # By convention, these are BedFile columns in a strict order.
  bed.cols <- c("chrom", "chromStart", "chromEnd", "name", "score",
                "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
                "blockSizes", "blockStarts")

  lim = 30000000
  row.num <- length(count.fields(bed.path))
  chunks.num <- row.num %/% lim
  skip.lens <- 0:chunks.num * lim

  for (skip.len in skip.lens) {

    bed <- fread(bed.path, skip = skip.len, nrow = lim, showProgress = FALSE)
    setnames(bed, names(bed), bed.cols[1:ncol(bed)])
    bed[, chromStart := chromStart + 1]

    # Just to set it apart from bedfile
    setnames(bed, c("chromStart", "chromEnd"), c("start", "end"))

    # Check for strand polarity
    has.strand <- "strand" %in% names(bed)
    if (has.strand) {
      strand.type <- bed[, unique(strand)]
      strand.sensitive <- ifelse(length(strand.type) == 1, FALSE, TRUE)
    } else {
      strand.sensitive <- FALSE
    }

    # Calculate unique length(s)
    len <- bed[, unique(end - start + 1)]

    dir.create(output.path, recursive = TRUE, showWarnings = FALSE)

    bed[, fwrite(.SD,  paste0(output.path, "/", chrom, ".csv",
                              if(compress) ".gz"),
                 append = TRUE),

       # .SDcols = c("start",
       #             if (length(len) > 1) "end",
       #             if (strand.sensitive) "strand"),
        by = chrom]

    # Save attribute
    att <- data.table(is.strand.sensitive = strand.sensitive)
    if (has.strand & !strand.sensitive) att$reported.strand <- strand.type
    if (length(len) == 1) att$case.length <- len

    fwrite(att, paste0(output.path, "/attributes.csv"), showProgress = FALSE)

  }
}
