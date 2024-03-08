#' Function retrieves the latest version information of the COSMIC database 
#' and the associated genome version by scraping data from the COSMIC website.
#'
#' @return A named vector containing the latest COSMIC version (`cosmic`) and 
#'         genome version (`genome`).
#'
#' @importFrom data.table fread
#' @importFrom stringi stri_extract_first_regex
#' 
#' @export
getCOSMIClatestVersion <- function() {
  
  scrap <- fread("https://cancer.sanger.ac.uk/census", header = FALSE,
                 sep = "", quote = "", showProgress = FALSE)[
                   V1 %like% "COSMIC v[0-9]+" & V1 %like% "GRCh"]
  cosmic.version <- scrap[, stri_extract_first_regex(V1, "v[0-9]+")]
  genome.version <- scrap[, stri_extract_first_regex(V1, "GRCh[0-9]+")]
  
  return(c(cosmic = cosmic.version, genome = genome.version))
}