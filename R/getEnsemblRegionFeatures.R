#' Get features of a given region.
#' 
#' GET overlap/region/:species/:region?feature=XXX?feature=XXXX?....
#' 
#' somatic_variation is from COSMIC
#'
#' @param species A species name or alias e.g. homo_sapiens, human, etc.
#' @param chromosome Ensembl chromosome name. Unlike UCSC, no prefix chr.
#' @param start Start position.
#' @param end End position.
#' @param features Region features. Options are "band", "gene", "transcript",
#'    "cds", "exon", "repeat", "simple", "misc", "variation",
#'    "somatic_variation", "structural_variation",
#'    "somatic_structural_variation", "constrained", "regulatory", "motif",
#'    "peak", "other_regulatory", "array_probe", and "mane".
#' @return A `data.table` of Ensembl features.
#'
#' @export
getEnsemblRegionFeatures <- function(species, chromosome, start, end,
                                     features) {
  
  feature.list <- c("band", "gene", "transcript", "cds", "exon", "repeat",
                    "simple", "misc", "variation", "somatic_variation",
                    "structural_variation", "somatic_structural_variation",
                    "constrained", "regulatory", "motif", "peak",
                    "other_regulatory", "array_probe", "mane")
  
  if (any(!features %in% feature.list)) {
    unknown.features <- features[!features %in% feature.list]
    stop("Unknown features: ", unknown.features)
  }
  
  server.limit <- 5e6 # 5Mb base
  if (end - start > server.limit) stop("Region is limited to 5Mb.")
  
  region <- paste0(chromosome, ":", start, "-", end)
  
  # Multiple features are accepted.
  features <- as.list(features)
  names(features) <- rep("feature", length(features))
  
  url <- paste0("http://rest.ensembl.org/overlap/region/", species, "/", region)
  url <- buildRESTurl(url, .list = features)
  
  # Server allow several output format e.g. bed, gff3, xml but fix format to
  # json only.
  h <- curl::new_handle()
  curl::handle_setheaders(h, "content-type" = "application/json")
  
  feature.dt <- getEnsemblData(url, h)
  setDT(feature.dt)
  
  return(feature.dt)
}