getCOSMIClatestVersion <- function() {
  
  scrap <- fread("https://cancer.sanger.ac.uk/census", header = FALSE,
                 sep = "", quote = "", showProgress = FALSE)[
                   V1 %like% "COSMIC v[0-9]+" & V1 %like% "GRCh"]
  cosmic.version <- scrap[, stri_extract_first_regex(V1, "v[0-9]+")]
  genome.version <- scrap[, stri_extract_first_regex(V1, "GRCh[0-9]+")]
  
  return(c(cosmic = cosmic.version, genome = genome.version))
}