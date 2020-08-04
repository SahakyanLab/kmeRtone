downloadCancerGeneCensus <- function(username, password, output) {

  id <- paste0(username, ":", password)

  raw.id <- charToRaw(id)

  auth.header <- paste("Basic", base64encode(raw.id))
  names(auth.header) <- "Authorization"

  download.file(

      url = paste0("https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/",
                   "cosmic/v91/cancer_gene_census.csv"),
      destfile = paste0(output , "cancer_gene_census.csv"),
      headers = auth.header,
      quiet = TRUE

  )

  url <- strsplit(readLines(paste0(output, "cancer_gene_census.csv"), warn=F),
                  '"')[[1]][4]

  download.file(url, destfile = paste0(output, "cancer_gene_census.csv"),
                quiet = TRUE)

}
