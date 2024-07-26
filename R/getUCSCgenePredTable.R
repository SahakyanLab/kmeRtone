#' Retrieve Gene Prediction Table from UCSC for a Given Genome
#'
#' This function retrieves the gene prediction table from the UCSC genome database
#' for a specified genome. It can fetch data from either the RefSeq or GENCODE databases.
#'
#' @param genome.name A string specifying the UCSC genome name for which the gene 
#'    prediction table is to be retrieved, e.g., 'hg38', 'mm39'.
#' @param db A string specifying the database used by UCSC to generate the table. 
#'    Options are 'refseq' or 'gencode'.
#'
#' @return A `data.table` containing the gene prediction table from the specified
#'    UCSC genome and database.
#'
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table setDT rbindlist
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