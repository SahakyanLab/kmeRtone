#' Select genomes for cross-species studies.
#'
#' The following filters are applied:
#' 1) assembly_level: "Complete Genome" or "Chromosome"
#' 2) genome_rep: "Full"
#' 3) Unique species_taxid (single representative species)
#' 4) refseq_category of "reference genome" is prioritised over
#'    "representative genome"
#'
#' @param organism.group Species group: archaea, bacteria, fungi, invertebrate,
#'    plant, protozoa, vertebrate_mammalian, vertebrate_other, or viral.
#' @param db Database record to use: refseq or genbank
#'
#' @return NCBI assembly summary with added column organism.group.
selectGenomesForCrossSpeciesStudy <- function(organism.group = "bacteria",
                                              db = "refseq") {

  asm <- getNCBIassemblySummary(organism.group = organism.group, db = db)

  # Select complete full genome only
  asm <- asm[assembly_level %in% c("complete genome", "chromosome") &
               genome_rep == "full"]

  # Sort to prioritize reference over representative genome
  asm[refseq_category == "reference genome", temp := 1]
  asm[refseq_category == "representative genome", temp := 2]
  setorder(asm, temp, na.last = TRUE)[, temp := NULL]

  # Only retain single species representative i.e. remove multiple strains
  asm <- asm[!duplicated(species_taxid)]

  # Assign column organism.group
  asm[, organism_group := organism.group]

  return(asm)
}