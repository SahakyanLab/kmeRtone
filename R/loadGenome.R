#' Build Genome object.
#'
#' The Genome object is capable of loading chromosome sequence on demand.
#'     UCSC Genomes are included in this kmeRtone package. Their specific
#'     chromosome sequence will be downloaded on demand once.
#'
#' @param genome.name A genome name. UCSC and NCBI genome is included with
#'    kmeRtone. Input their name e.g. hg19 or GRCh37.
#' @param fasta.style FASTA version: "UCSC" or "NCBI".
#' @param ncbi.db NCBI database: "refseq" or "genbank".
#' @param ncbi.asm NCBI assembly table.
#' @param mask Genome mask: "none", "soft", or "hard". Default is "none".
#' @param fasta.path Path to a directory of user-provided genome FASTA files or 
#'  the destination to save the NCBI/UCSC downloaded reference genome files.
#' @param use.UCSC.name For NCBI Genome, use UCSC-style chromosome name? Default
#'    is FALSE.
#' @param load.limit Maximum chromosome sequences loaded. Default is 1.
#' @return A `UCSC_Genome` or `NCBI_Genome` object.
#'
#' @export
loadGenome <- function(genome.name, fasta.style, mask="none", fasta.path,
                       ncbi.db, ncbi.asm, use.UCSC.name=FALSE, load.limit=1) {

  if (fasta.style == "UCSC") {
    genome <- UCSC_Genome$new(genome.name, fasta.path, mask, load.limit)
  } else if (fasta.style == "NCBI") {
    genome <- NCBI_Genome$new(genome.name, ncbi.db, fasta.path, ncbi.asm, mask,
                              use.UCSC.name, load.limit)

  }

  return(genome)
}
