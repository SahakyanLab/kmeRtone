#' Function downloads chromosome-separated fasta genome sequences from the UCSC
#' database. Users can specify a genome name, an output folder, and a specific
#' chromosome or chromosomes. There's an option to choose the download method as well.
#'
#' @param genome.name Genome name (e.g., hg19, hg38, mm19).
#' @param output.path Output folder for the downloaded sequences.
#' @param chr.name Specific chromosome to download; defaults to all if unspecified.
#' @param method Download method for the `download.file` function.
#'
#' @return An output folder containing chromosome-separated fasta files.
#'
#' @importFrom utils download.file
#' @importFrom data.table fread fwrite
#' @importFrom stringi stri_extract_all_regex
#' 
#' @export
downloadUCSCgenome <- function(genome.name, output.path, chr.name,
                               method = "curl") {

  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  fasta.pattern <- "\\.(fa|fna|fasta)($|\\.gz$)"

  root.url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/", genome.name,
                     "/chromosomes/")
  root.page <- fread(root.url, sep = "", showProgress = FALSE)

  chr.files <-
    stri_extract_all_regex(root.page[[1]],
                           "chr([0-9]{1,2}|[XYM]|[IVX]{1,4})\\.fa\\.gz",
                           omit_no_match = TRUE) |> unlist() |> unique()

  # Save all available chromosomes in info.csv
  if (!file.exists(paste0(output.path, "/info.csv")))
    fwrite(data.table(chromosome = sub(fasta.pattern, "", chr.files) |>
                                   rep(each = 2),
                      strand = c("+", "-")),
           paste0(output.path, "/info.csv"))

  # Check for chromosome existence
  if (missing(chr.name)) {

    chr.name <- sub(fasta.pattern, "", chr.files)

  } else {

    miss.files <- chr.name[!chr.name %in% sub(fasta.pattern, "", chr.files)]
    if (length(miss.files) == length(chr.name)) {
      stop("No chromosome name: ", miss.files)
    } else if (length(miss.files) > 0) {
      warning("No chromosome name: ", miss.files)
    }

  }

  chr.files <- chr.files[sub(fasta.pattern, "", chr.files) %in% chr.name]

  for (chr.file in chr.files) {

    chr.path <- paste0(output.path, "/", chr.file)
    chr.url <- paste0(root.url, "/", chr.file)

    cat("\nDownloading", chr.file, "from UCSC...\n")
    download.file(chr.url, chr.path, method = method)

  }

}
