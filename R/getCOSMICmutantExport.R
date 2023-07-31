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