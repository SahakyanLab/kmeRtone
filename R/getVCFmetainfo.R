#' Read VCF metainfo file using tabix.
#'
#' Require tabix in PATH
#' VCF manual is referred from https://samtools.github.io/hts-specs/VCFv4.3.pdf
#'
#' @param vcf.file A path to a local or remote tabix-indexed VCF file.
#'
#' @return VCF metainfo.
#'
#' @export
getVCFmetainfo <- function(vcf.file) {
  system3("tabix", args = c("-DH", vcf.file), stdout = TRUE)
}