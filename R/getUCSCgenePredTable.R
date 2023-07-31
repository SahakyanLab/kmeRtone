#' Study k-mer composition across species.
#'
#' @param UCSC.genome.name UCSC genome name e.g. hg38 and mm39.
#' @param db Database used by UCSC to generate the table: "refseq" or "gencode"
#'
#' @return genepred `data.table`.
#'
#' @export
getUCSCgenePredTable <- function(genome.name, db) {
  
  table.name <- if(db == "refseq") "ncbiRefSeq"
                else if (db == "gencode") "knownGene"
  rest.url <- "https://api.genome.ucsc.edu/getData/track"
  genepred.url <- buildRESTurl(url = rest.url, genome = genome.name,
                               track = table.name)
  response <- curl::curl_fetch_memory(genepred.url)
  genepred <- jsonlite::fromJSON(rawToChar(response$content))[[table.name]]
  
  if (db == "refseq") genepred <- rbindlist(genepred)
  else setDT(genepred)

  return(genepred)
}