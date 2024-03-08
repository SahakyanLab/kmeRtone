#' Get gnomAD VCF file using tabix.
#' 
#' Function retrieves variant data from gnomAD VCF files using tabix for a
#' specified set of genomic regions. It allows users to select the gnomAD
#' version and server location (Google, Amazon, or Microsoft) for fetching the
#' data.
#'
#' @param chr.names Chromosome names.
#' @param starts Start positions.
#' @param ends End positions.
#' @param INFO.filter Parse only filtered INFO ID. Default is to parse all IDs.
#' @param version The gnomAD version. Default to latest version 3.1.2.
#' @param server Server locations: "google", "amazon", or "microsoft". Default
#'    is random.
#'
#' @return A data.table of VCF.
#' 
#' @importFrom data.table rbindlist
#'
#' @export
getGnomADvariants <- function(chr.names, starts, ends, INFO.filter=NULL,
                              version="3.1.2", server="random") {

  # The remote VCF files are separated by chromosome XXX
  # In the future version, the url might break.
  google.vcf.url <-
    paste0("https://storage.googleapis.com/gcp-public-data--gnomad/release/",
           version, "/vcf/genomes/gnomad.genomes.v", version,
           ".sites.XXX.vcf.bgz")

  aws.vcf.url <- paste0("https://gnomad-public-us-east-1.s3.amazonaws.com/",
                        "release/", version, "/vcf/genomes/gnomad.genomes.v",
                        version, ".sites.XXX.vcf.bgz")

  # It gives error: Illegal seek. Something is wrong with the remote file.
  microsoft.vcf.url <-
    paste0("https://datasetgnomad.blob.core.windows.net/dataset/release/",
           version, "/vcf/genomes/gnomad.genomes.v", version,
           ".sites.XXX.vcf.bgz")

  vcf.dt <- lapply(unique(chr.names), function(chr.name) {

    if (server == "google") {
      vcf.url <- google.vcf.url
    } else if (server == "amazon") {
      vcf.url <- aws.vcf.url
    } else if (server == "microsoft") {
      vcf.url <- microsoft.vcf.url
    } else if (server == "random") {
      vcf.url <- sample(c(google.vcf.url, aws.vcf.url), 1)
    }
    vcf.file <- vcf.url <- sub("XXX", chr.name, vcf.url)

    idx <- which(chr.names == chr.name)
    vcf.dt <- readVCF(vcf.url, chr.name, starts[idx], ends[idx], INFO.filter)

  }) |> rbindlist()

  metainfo <- getVCFmetainfo(vcf.file)
  class(vcf.dt) <- c("VCF", class(vcf.dt))
  setattr(vcf.dt, "metainfo", metainfo)

  return(vcf.dt)
}