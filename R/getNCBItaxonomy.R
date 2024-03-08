#' Get all taxonomy classification from NCBI
#'
#' NCBI provides taxonomy for all life domains in rankedlineage.dmp file.
#' This function download and store the file in kmeRtone_data path.
#'
#' @return A data.table of taxonomy classification.
#' 
#' @importFrom tools md5sum
#' @importFrom data.table fread set setnames
#' @importFrom utils head tail untar
#' 
#' @export
getNCBItaxonomy <- function() {

  taxonomy.dir <- paste0(path.expand(kmeRtone.data.path), "/taxonomy/")
  file.url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
  file.name <- paste0(taxonomy.dir, "new_taxdump.tar.gz")
  md5.url <- paste0(file.url, ".md5")

  # Check if new_taxdump.tar.gz is already exist and its md5sum
  if (file.exists(file.name)) {
    md5.ftp <- fread(md5.url, header = FALSE, showProgress = FALSE)[[1]]
    md5.local <- tools::md5sum(file.name)
  }

  # Download new_taxdump.tar.gz
  if (!file.exists(file.name) || md5.ftp != md5.local) {
    persistentDownload(file.url, file.name)
  }

  # Extract only rankedlineage.dmp
  if (!file.exists(tax.file <- paste0(taxonomy.dir, "rankedlineage.dmp"))) {
    untar(file.name, files = "rankedlineage.dmp", exdir = taxonomy.dir)
  }

  tax <- fread(tax.file, showProgress = FALSE)

  # Rename columns
  setnames(tax, paste0("V", seq(1, 20, 2)),
           c("tax_id", "tax_name", "species", "genus", "family", "order",
             "class", "phylum", "kingdom", "superkingdom"))

  # Remove | separator
  set(tax, j = paste0("V", seq(2, 20, 2)), value = NULL)

  return(tax)
}
