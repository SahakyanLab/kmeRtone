#' Function downloads the latest Cosmic Mutant Export data from the COSMIC
#' database. It requires the user to be registered with COSMIC for non-commercial
#' use. The function constructs the URL for the latest mutant export file, 
#' authenticates the URL, and then downloads the data.
#'
#' @param email Email registered with COSMIC for accessing data.
#' @param password Password for the COSMIC account.
#'
#' @return A `data.table` containing the Cosmic Mutant Export data.
#'
#' @importFrom data.table fread
#' 
#' @export
getCOSMICmutantExport <- function(email, password) {
  
  vers <- getCOSMIClatestVersion()
  cosmic.version <- vers["cosmic"]
  genome.version <- vers["genome"]
  
  mutants.url <- paste0("https://cancer.sanger.ac.uk/cosmic/file_download/",
                    genome.version, "/cosmic/", cosmic.version,
                    "/CosmicMutantExport.tsv.gz")
  
  mutants.url <- getCOSMICauthURL(email, password, mutants.url)
  
  mutants <- fread(mutants.url, showProgress = FALSE)
  
  return(mutants)
}