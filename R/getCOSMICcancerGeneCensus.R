#' Get Cancer Gene Census (CGC) from COSMIC database.
#' 
#' To access the data for non-commercial usage, you must register with the
#'    COSMIC. This function fetch the latest CGC.
#'
#' @param email Email registered with COSMIC.
#' @param password Password associated with the registered email.
#'
#' @return A `data.table` containing the Cancer Gene Census data.
#'
#' @importFrom data.table fread
#' @importFrom utils download.file
#' 
#' @export
getCOSMICcancerGeneCensus <- function(email, password) {
  
  vers <- getCOSMIClatestVersion()
  cosmic.version <- vers["cosmic"]
  genome.version <- vers["genome"]
  
  cgc.url <- paste0("https://cancer.sanger.ac.uk/cosmic/file_download/",
                    genome.version, "/cosmic/", cosmic.version,
                    "/cancer_gene_census.csv")
  
  cgc.url <- getCOSMICauthURL(email, password, cgc.url)
  
  cgc <- fread(cgc.url, showProgress = FALSE)
  
  return(cgc)
}