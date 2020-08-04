downloadHumanGenomeAnnotation <- function(assembly, track.table,
                                          output.dir=NULL) {
  # The flags reflect the Table Browser page
  # (http://genome.ucsc.edu/cgi-bin/hgTables) at UCSC Genome Browser.
  # This function only applicable to download human genome annotation.
  # Description of flags are as below:
  #   assembly       : Specifies which version of the organism's genome sequence
  #                    to use. Option: {"hg19", "hg38"}
  #   track.table    : Selects the annotation track data to work with. For
  #                    further information, refer to the URL above.
  #                               Track              Table
  #                    Option: NCBI RefSeq  {"ncbiRefSeq", "refGene"}
  #                            UCSC Genes   {"knownGene"}
  #                    For other choices, please refer to the "table" selectio
  #                    in the URL above.
  #   output.folder  : A folder path to save the genome annotation. Default i
  #                    current working directory.

  if (!assembly %in% c("hg19", "hg38")) {
    stop(paste("Assembly name", assembly, "is not recognized. Only hg19 and",
               "hg38 are supported."))
  }

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  ftp.url = paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/", assembly,
                   "/database/", track.table, ".txt.gz")
  save.path <- paste0(output.dir, "/", track.table, ".txt.gz")

  # Download the genome annotation
  download.file(ftp.url, save.path, quiet = TRUE)

  cat("\n-- Genome annotation table is saved at", save.path, "\n")

}
