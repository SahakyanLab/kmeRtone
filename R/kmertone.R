#' kmeRtone all-in-one user interface
#'
#' @description
#' This function serves as an all-in-one interface for various genomic data analyses
#' leveraging k-mer based techniques.
#'
#' @param case.coor.path Path to a folder containing chromosome-separated coordinate
#'        files or bedfiles. Assumed replicates for subfolder or bedfiles.
#' @param genome.name Name of the genome (e.g., "hg19", "hg38"). Default is "unknown".
#' @param strand.sensitive Logical value indicating whether strand polarity matters. 
#'        Default is TRUE.
#' @param k Length of k-mer to be investigated. Recommended values are 7 or 8.
#' @param ctrl.rel.pos A vector of two integers specifying the relative range positions 
#'        of control regions.
#' @param case.pattern Regular expression pattern for identifying case regions.
#'        Default is NULL.
#' @param output.dir Directory path for saving output files. Default is "output/".
#' @param case Optional pre-built Coordinate object.
#' @param genome Optional pre-built Genome object.
#' @param control Optional pre-built control Coordinate object.
#' @param control.path Path for pre-built control Coordinate object.
#' @param genome.path Path to a directory of user-provided genome FASTA files.
#' @param rm.case.kmer.overlaps Logical indicating whether to remove overlapping 
#'        k-mers in case regions. Default is FALSE.
#' @param single.case.len Integer indicating uniform length of case regions.
#' @param merge.replicates Logical indicating whether to merge replicates. 
#'        Default is TRUE.
#' @param kmer.table Pre-calculated k-mer score table.
#' @param module Selected kmeRtone module to run. Possible values include "score", 
#'        "explore", "tune", among others.
#' @param rm.dup Logical indicating whether to remove duplicate coordinates.
#'        Default is TRUE.
#' @param case.coor.1st.idx Integer specifying indexing format for case coordinates.
#' @param ctrl.coor.1st.idx Integer specifying indexing format for control coordinates.
#' @param coor.load.limit Maximum number of coordinates to load. Default is 1.
#' @param genome.load.limit Maximum number of genome sequences to load. Default is 1.
#' @param genome.fasta.style String specifying the style of the genome FASTA. 
#'        Possible values are "UCSC", "NCBI". Default is "UCSC".
#' @param genome.ncbi.db String specifying the NCBI database to use. Possible values
#'        are "refseq", "genbank". Default is "refseq".
#' @param use.UCSC.chr.name Logical indicating whether to use UCSC chromosome names.
#' @param verbose Logical indicating whether to display progress messages.
#'        Default is TRUE.
#' @param kmer.cutoff Cutoff percentage for k-mer selection in case studies.
#'        Default is 5.
#' @param selected.extremophiles Vector of selected extremophile species for study.
#' @param other.extremophiles Vector of other extremophile species for control.
#' @param cosmic.username COSMIC username for accessing the cancer gene census.
#' @param cosmic.password COSMIC password for accessing the cancer gene census.
#' @param tumour.type.regex Regular expression pattern for filtering tumour types.
#' @param tumour.type.exact Exact tumour type to be included in the cancer gene census.
#' @param cell.type Cell type to be included in the cancer gene census. Default is
#'        "somatic".
#' @param genic.elements.counts.dt Data table of susceptible k-mer counts in genic
#'        elements.
#' @param population.size Size of the population for cross-population studies. 
#'        Default is 1 million.
#' @param selected.genes Selected genes for mutation in cross-population studies.
#' @param add.to.existing.population Logical indicating whether to add to the existing
#'        simulated population. Default is FALSE.
#' @param population.snv.dt Data table of single nucleotide variants used in 
#'        population simulations.
#' @param pop.plot Logical indicating whether to plot the outcome of the cross-population
#'        study. Default is TRUE.
#' @param pop.loop.chr Logical indicating whether to loop based on chromosome name 
#'        in cross-population studies. Default is FALSE.
#'
#' @return Depends on the selected module.
#'
#' @export
kmeRtone <- function(case.coor.path, genome.name, strand.sensitive, k,
                     ctrl.rel.pos=c(80, 500), case.pattern,
                     output.dir="output/", case, genome, control, control.path,
                     genome.path, rm.case.kmer.overlaps, single.case.len,
                     merge.replicates, kmer.table, module="score", rm.dup=TRUE,
                     case.coor.1st.idx=1, ctrl.coor.1st.idx=1,
                     coor.load.limit=1, genome.load.limit=1,
                     genome.fasta.style="UCSC", genome.ncbi.db="refseq",
                     use.UCSC.chr.name=FALSE,
                     verbose=TRUE, kmer.cutoff=5, selected.extremophiles,
                     other.extremophiles, cosmic.username, cosmic.password,
                     tumour.type.regex=NULL, tumour.type.exact=NULL,
                     cell.type="somatic", genic.elements.counts.dt,
                     population.size=1e6, selected.genes,
                     add.to.existing.population=FALSE, population.snv.dt=NULL,
                     pop.plot=TRUE, pop.loop.chr=FALSE) {

  # Argument checking ----------------------------------------------------------
  message("\n")

  # Coordinate
  if (missing(case.coor.path) & missing(case) &
      module %in% c("score", "explore")) {
    stop("Please provide case coordinate information.")
  }
  if (!missing(case)) stopifnot("Coordinate" %in% class(case))

  # Genome
  if (missing(genome.name) & missing(genome) & missing(genome.path) &
      module %in% c("score", "explore", "study_genic_elements")) {
    stop("Please provide genome information.")
  }
  if (!missing(genome)) stopifnot("Genome" %in% class(genome))

  # Operation
  if ((missing(strand.sensitive) || !is.logical(strand.sensitive)) &
      module %in% c("score", "explore")) {
    stop("Please indicate whether strand.sensitive is TRUE or FALSE.")
  }
  if ((missing(k) || !is.numeric(k)) &
      ! module %in% c("explore", "study_cancer_genes")) {
    stop("Please indicate length of k.")
  }
  if ((missing(rm.case.kmer.overlaps) || !is.logical(rm.case.kmer.overlaps)) &
      module == "score") {
    stop("Please indicate whether to rm.case.kmer.overlaps or not. TRUE/FALSE")
  }
  if (missing(case) & (missing(merge.replicates) ||
                       !is.logical(merge.replicates)) &
      module %in% c("score")) {
    if (length(list.files(case.coor.path,
                          "^chr.+\\.(csv|txt|tsv)($|\\.gz$)")) == 0) {
      stop("Please indicate whether to merge.replicates or not. TRUE/FALSE")
    } else {
      merge.replicates <- FALSE
    }
  }
  if (missing(kmer.table) & module == "score") {
    kmer.table <- NULL
  } else if (missing(kmer.table) & module %in% c("study_across_species",
      "study_genic_elements")) {
        stop("Please input kmer.table for table of k-mer scores.")
      }
  if (!missing(case.pattern) && !is.null(case.pattern) &&
      any(stri_detect_regex(case.pattern, "[ACGT]", negate = TRUE))) {
    stop("Only base A, C, G, and T are supported.")
  }

  # Pre-populate missing variable
  if (missing(case)) case <- NULL
  if (missing(case.coor.path)) case.coor.path <- NULL
  if (missing(genome)) genome <- NULL
  if (missing(genome.name)) genome.name <- "unknown"
  if (missing(genome.path)) genome.path <- NULL
  if (missing(control)) control <- NULL
  if (missing(control.path)) control.path <- NULL
  if (missing(single.case.len)) {
    if (missing(case)) {
      message(paste("Argument single.case.len is not specified. Case length is assumed",
          " to be varied.\n"))
      single.case.len <- NULL
    } else {
      single.case.len <- case$single_len
    }
  }
  if (missing(case.pattern) & module != "study_cancer_genes") {
    message(paste("Case pattern is not specified. Case pattern is assumed to be any",
        "pattern i.e. it is set to NULL.\n"))
    case.pattern <- NULL
  }

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  # MAIN -----------------------------------------------------------------------

  if ("score" %in% module) {

    kmer.table <- SCORE(
      case.coor.path, genome.name, strand.sensitive, k, ctrl.rel.pos,
      case.pattern, output.dir, case, genome, control, control.path,
      genome.path, rm.case.kmer.overlaps, single.case.len, merge.replicates,
      rm.dup, case.coor.1st.idx, ctrl.coor.1st.idx, coor.load.limit,
      genome.load.limit, genome.fasta.style, genome.ncbi.db, use.UCSC.chr.name,
      verbose)

    if (length(module) == 1) return(kmer.table)

  }

  if ("explore" %in% module) {

    EXPLORE(
      case.coor.path, genome.name, strand.sensitive, k, case.pattern,
      output.dir, case, genome, control, genome.path, single.case.len, rm.dup,
      case.coor.1st.idx, coor.load.limit, genome.load.limit, genome.fasta.style,
      genome.ncbi.db, use.UCSC.chr.name, verbose)
  }

  if ("tune" %in% module) TUNE()

  if ("study_across_species" %in% module) {

    STUDY_ACROSS_SPECIES(
      kmer.table, kmer.cutoff, k, case.pattern, selected.extremophiles,
      other.extremophiles, paste0(output.dir, "/study_across_species"))

  }

  if ("study_genic_elements" %in% module) {
    STUDY_GENIC_ELEMENTS(
      kmer.table, kmer.cutoff, k, genome.name, case.pattern,
      db = "refseq", paste0(output.dir, "/study_genic_elements/"))
  }

  if ("study_cancer_genes" %in% module) {
    STUDY_CANCER_GENES(cosmic.username, cosmic.password, tumour.type.regex,
      tumour.type.exact, cell.type, genic.elements.counts.dt,
      paste0(output.dir, "/study_cancer_genes/"))
  }

  if ("study_across_populations" %in% module) {
    STUDY_ACROSS_POPULATIONS(kmer.table, kmer.cutoff, genome.name, k,
      db = "refseq", case.pattern, population.size, selected.genes,
      add.to.existing.population = add.to.existing.population,
      paste0(output.dir, "/study_across_populations/"),
      population.snv.dt, pop.loop.chr, pop.plot)
  }

  if ("study_G4_susceptibility" %in% module) STUDY_G4_SUSCEPTIBILITY()

}