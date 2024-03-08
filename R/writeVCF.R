#' Write VCF file and compress using bgzip.
#' 
#' Require bgzip in PATH
#' VCF manual is referred from https://samtools.github.io/hts-specs/VCFv4.3.pdf
#'
#' @param vcf A VCF object.
#' @param output.vcf.gz Output filename including vcf.gz extension.
#' @param append To append or not? Default is FALSE.
#' @param tabix To tabix or not? Default is FALSE.
#'
#' @importFrom data.table set setorder fwrite setalloccol
#' @importFrom stringi stri_paste stri_sub_replace stri_split_fixed stri_extract_first_regex stri_replace_first_fixed stri_replace_all_regex
#'
#' @export
writeVCF <- function(vcf, output.vcf.gz, append=FALSE, tabix=FALSE) {
  
  if (!file.exists(output.vcf.gz) & append) append <- FALSE

  tmp.vcf <- tempfile(fileext = ".vcf")
  
  metainfo <- attr(vcf, "metainfo")
  
  if (!append) writeLines(attr(vcf, "metainfo"), tmp.vcf)
  
  info.ids <- metainfo[grepl("##INFO", metainfo)] |>
    stri_extract_first_regex("ID=.+?,") |> stri_replace_all_regex("^ID=|,$", "")
  
  # Check if there is no more allocation for data.table column
  if (length(vcf) >= truelength(vcf)) setalloccol(vcf, truelength(vcf) + 1)
  
  set(vcf, j = "INFO", value = "")
  
  for (id in info.ids) {
    
    info <- vcf[[paste0("INFO.", id)]]
    
    if (is.list(info)) {
      info <- lapply(info, function(info) {
        info[is.na(info)] <- "."
        info <- paste(info, collapse = ",")
        info
      })
    }
    
    if (is.logical(info)) {
      idx <- which(info)
      info[idx] <- id
      info[-idx] <- ""
    } else {
      info <- stri_paste(id, "=", info)
    }
    
    set(vcf, j = "INFO", value = stri_paste(vcf$INFO, ";", info))
    
  }
  
  # Remove ; at front.
  set(vcf, j = "INFO",
      value = stri_sub_replace(vcf$INFO, from = 1, to = 1, replacement = ""))
  
  header <- (metainfo[length(metainfo)] |> stri_replace_first_fixed("#", "") |>
    stri_split_fixed("\t"))[[1]]
  
  setorder(vcf, CHROM, POS)
  
  vcf[, fwrite(.SD, tmp.vcf, sep = "\t", append = TRUE, compress = "none",
               col.names = FALSE, quote = FALSE, showProgress = FALSE),
      .SDcols = header]
  
  set(vcf, j = "INFO", value = NULL)
  
  system(paste0("bgzip -c ", tmp.vcf, if (append) " >> " else " > ",
                output.vcf.gz))

  if (tabix) system(paste0("tabix -fp vcf ", output.vcf.gz))
  
  invisible(NULL)
}