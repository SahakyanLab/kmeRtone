#' Get features of given variant IDs.
#' 
#' POST variation/:species
#'
#' @param species Species name or alias e.g. homo_sapiens, human, etc.
#' @param variant.ids A vector of variant IDs. e.g. rs56116432, COSM476, etc.
#' @param include.genotypes Include genotypes? Default is FALSE.
#' @param include.phenotypes Include phenotypes? Default is FALSE.
#' @param include.allele.frequencies Include allele frequencies? Default is
#'    FALSE.
#' @param include.genotype.frequencies Include genotype frequencies? Default is
#'    FALSE.
#' @return A variant-named list of list of variation features.
#'
getEnsemblVariantFeatures_serial <- function(species, variant.ids,
                                      include.genotypes=FALSE,
                                      include.phenotypes=FALSE,
                                      include.allele.frequencies=FALSE,
                                      include.genotype.frequencies=FALSE) {
  
  # Build REST URL
  url <- paste0("https://rest.ensembl.org/variation/", species)
  includes <- c(include.genotypes, include.phenotypes,
                include.allele.frequencies, include.genotype.frequencies) |>
    as.integer() |> as.list()
  names(includes) <- c("genotypes", "phenotypes", "pops",
                       "population_genotypes")
  url <- buildRESTurl(url, .list = includes)
  
  h <- curl::new_handle()
  curl::handle_setheaders(h, "content-type" = "application/json")
  
  # Server limit to 200 variant IDs. For more than 200 IDs, loop.
  feature.lists <- lapply(seq(1, length(variant.ids), 200), function(i) {
    
    variant.ids.json <- paste0('{"ids" : ["',
                               paste(variant.ids[i:min(length(variant.ids),
                                                       i+200-1)],
                                     collapse = '", "'),
                               '" ]}')
    
    curl::handle_setopt(h, customrequest = "POST",
                        postfields = variant.ids.json)
    
    feature.lists <- getEnsemblData(url, h)
    
    feature.lists
  }) |> unlist(recursive = FALSE)
  
  return(feature.lists)
}