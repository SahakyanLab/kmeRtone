#' Function downloads genome fasta files from the NCBI FTP database. Users can
#' provide either organism names or an assembly summary data table. 
#' 
#' Supports options for splitting multi-header fasta files and overwriting existing files.
#'
#' @param species Species names.
#' @param asm NCBI assembly summary data.table
#' @param db Database record to use: refseq or genbank
#' @param output.dir Output directory path. Default is current directory.
#' @param split.fasta NCBI fasta files are multi-header. Split them? Default is
#'    FALSE.
#' @param overwrite Overwrite any existed genome file? Default is FALSE to skip
#'    the download.
#'
#' @return Genome fasta file(s) named according to the FTP database convention.
#'
#' @importFrom data.table data.table
#' @importFrom stringi stri_split_regex
#' 
#' @export
downloadNCBIGenomes <- function(asm, species, db, output.dir="./",
                                split.fasta=FALSE, overwrite=FALSE) {

  if (!missing(species) & missing(asm)) {
    asm <- getNCBIassemblySummary(organism.group = "all", db = db)
    asm <- asm[organism_name %in% species]
    asm <- selectRepresentativeFromASM(asm)
    if (nrow(asm) == 0)
      stop("Species name(s) not found in NCBI assembly summary")
  }

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  if (!overwrite) {
    asm <- asm[!file.exists(paste0(output.dir, basename(ftp_path),
                                   "_genomic.fna.gz"))]
    if (nrow(asm) == 0) return(invisible(NULL))
  }

  asm[, {

    # Try to guess link if na
    if (ftp_path == "na") {
      dir.db <- stri_split_regex(assembly_accession, "_|\\.")[[1]][1:2]
      dir.subs <- substring(dir.db[2], seq(1, nchar(dir.db[2]), 3),
                            seq(3, nchar(dir.db[2]), 3))
      dir.ftp <- paste(c(dir.db[1], dir.subs), collapse = "/")
      ftp_path <- paste0(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/", dir.ftp, "/",
        assembly_accession, "_", gsub("[^A-Za-z0-9]", "", asm_name))
    }

    # FASTA file
    file.name <- paste0(basename(ftp_path), "_genomic.fna.gz")
    url.path <- paste0(ftp_path, "/", file.name)
    file.path <- paste0(output.dir, "/", file.name)

    persistentDownload(file.url = url.path, output.name = file.path)

    if (split.fasta) {

      splitFASTA(file.path, output.dir = paste0(output.dir, "/",
                                                basename(ftp_path)))
      file.remove(file.path)
    }

    # Assembly report for accession ID details e.g. chromosome, plasmid, etc.
    file.name <- paste0(basename(ftp_path), "_assembly_report.txt")
    url.path <- paste0(ftp_path, "/", file.name)
    file.path <- paste0(output.dir, "/", file.name)

    persistentDownload(file.url = url.path, output.name = file.path)

  }, by = seq_len(nrow(asm))]

  return(invisible(NULL))
}
