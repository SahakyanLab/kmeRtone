#' Get NCBI assembly summary.
#'
#' @param db Database record to use: refseq or genbank
#' @param taxonomic.group Options are archaea, bacteria, fungi, invertebrate,
#'    plant, protozoa, vertebrate_mammalian, vertebrate_other, viral, or all.
#'
#' @return Assembly summary data.table.
#'
#' @export
getNCBIassemblySummary <- function(organism.group, db = "refseq") {

  group.selections <- c("archaea", "bacteria", "fungi", "invertebrate", "plant",
                        "protozoa", "vertebrate_mammalian", "vertebrate_other",
                        "viral")

  if (organism.group %in% group.selections) {
    asm.url <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/", db, "/",
                      organism.group, "/assembly_summary.txt")
  } else if (organism.group == "all") {
    asm.url <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/", db,
                      "/assembly_summary_", db, ".txt")
  } else {
    stop("organism.group should be either ",
         paste(group.selections, collapse = ", "), ", or all")
  }

  asm.tmp <- tempfile("asm", fileext = ".txt")
  persistentDownload(asm.url, asm.tmp)
  asm <- fread(asm.tmp, quote = "", showProgress = FALSE)

  # The labels used in asm table are not consistent in their letter case.
  # We use all small case for consistency.
  # Change to factor to sort
  set(asm, j = "refseq_category",
      value = stri_trans_tolower(asm$refseq_category))
  set(asm, j = "assembly_level",
      value = stri_trans_tolower(asm$assembly_level))
  set(asm, j = "genome_rep",
      value = stri_trans_tolower(asm$genome_rep))

  setnames(asm, colnames(asm)[1],
           stri_trim_left(sub("^#", "", colnames(asm)[1])))

  return(asm)
}