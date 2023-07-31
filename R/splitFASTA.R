#' Split a FASTA file by header.
#'
#' The first non-space word in the header will be used as the file name.
#'
#' data.table::fread is not built for reading in chunks. The developers expect
#' skip and nrow arguments to be in a small number. Large number slows the
#' reading a bit.
#'
#' @param fasta.file A path to a FASTA file.
#'
#' @return A splitted fasta files by its headers.
#'
#' @export
splitFASTA <- function(fasta.file, output.dir="./") {

  # Helper function for repetitive resolving filename
  getFilename <- function(header) {
    stri_extract_first_regex(header, "\\S+") |>
      stri_replace_first_fixed(">", "") |> paste0(".fa.gz")
  }

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  # To avoid loading all lines to RAM, read by chunk.
  # 12000000 row is ~2GB data.table for common width of fasta file.
  max.row <- 12000000

  skip.len <- 0
  chunk.num <- 0
  fasta.file.conn <- file(fasta.file, "rt")
  while (TRUE) {

    fasta <- scan(fasta.file, what = "character", sep = "\n", nlines = max.row,
                  quiet = TRUE, quote = "") |> as.data.table()

    # scan is faster. data.table developer did not expect skip to be large.
    # fasta <- fread(fasta.file, sep = "", showProgress = FALSE, header = FALSE,
    #                skip = skip.len, nrow = max.row)

    # Assign group number for every group of sequences
    fasta[, group := stri_startswith_fixed(V1, ">")][, group := cumsum(group)]

    # Get filename from the last header in the current chunk. Group 0 is has no
    # header and from previous chunk.
    if (fasta[1, group > 0])
      file.name <- fasta[group == max(group), getFilename(V1[1])]

    fasta[, {

      if (group > 0) file.name <- getFilename(V1[1])

      fwrite(.SD, paste0(output.dir, "/", file.name), showProgress = FALSE,
             col.names = FALSE, append = group == 0)

    }, by = group]

    chunk.num <- chunk.num + 1
    skip.len <- chunk.num * max.row
    if (nrow(fasta) < max.row) break
  }
  close(fasta.file.conn)

  return(invisible(NULL))
}
