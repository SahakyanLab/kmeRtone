#' Select the best representative species from the NCBI assembly summary.
#'
#' sort.idx is a weight to sort where heavier will be preffered. Any tie weight
#'    will be further sorted by organism_name. Only the top unique species_taxid
#'    will be retained in the final assembly summary.
#'
#' @param asm NCBI assembly summary.
#' @return Trimmed NCBI assembly summary.
#' @export
selectRepresentativeFromASM <- function(asm) {

  # The labels used in asm table are not consistent in their letter case.
  # We use all small case for consistency.
  refseq.categories <- c("reference genome", "representative genome", "na")
  assembly.levels <- c("complete genome", "chromosome", "scaffold", "contig",
    "na")
  genome.reps <- c("full", "partial", "na")

  # Change to factor to sort
  changes <- NULL
  set(asm, j = "refseq_category",
      value = factor(stri_trans_tolower(asm$refseq_category),
                     levels = refseq.categories))
  set(asm, j = "assembly_level",
      value = factor(stri_trans_tolower(asm$assembly_level),
                     levels = assembly.levels))
  set(asm, j = "genome_rep",
      value = factor(stri_trans_tolower(asm$genome_rep),
                     levels = genome.reps))

  # Sort in the following order to select the top as the representative
  setorder(asm, refseq_category, assembly_level, genome_rep)

  # Only retain single species representative i.e. remove multiple strains
  asm <- asm[!duplicated(species_taxid)]

  return(asm)
}