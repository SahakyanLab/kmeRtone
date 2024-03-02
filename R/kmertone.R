#' kmeRtone all-in-one user interface
#'
#' @field case.coor.path A path to a folder containing either (1) chromosome-
#'      separated coordinate files (assume replicates for subfolder) or (2)
#'      bedfile. (assume replicates for bedfiles).
#' @field genome.name Name of genome. kmeRtone include UCSC genome e.g. "hg19"
#'      and "hg38". Default is "unknown".
#' @field strand.sensitive Does strand polarity matters? Default is TRUE. In
#'      sensitive strand, case k-mers are extracted from the sense strand only.
#'      In contrast, in insensitive strand, the case k-mers are extracted from
#'      both strands.
#' @field k Length of k-mer to be investigated. Recommended is 7 or 8.
#' @field ctrl.rel.pos A relative range position of control regions. For example
#'      c(80,500) means control regions are 80-500 bases away from the the case
#'      site (upstream and downstream).
#' @field case.pattern Case pattern(s). Default is NULL for no pattern. Only
#'      applicable to module other than study_cancer_genes.
#' @field output.dir An output directory name. Default is output/.
#' @field case Optional pre-built Coordinate object. Option for
#'      strand.sensitive and single.case.len will be ignored.
#' @field genome Optional pre-built Genome object. Option genome.name and
#'      genome.path will be ignored.
#' @field control Optional pre-built control Coordinate object. Option
#'      ctrl.rel.pos will be ignored.
#' @field control.path A path for pre-built control Coordinate object. Option
#'      ctrl.rel.pos will be ignored.
#' @field genome.path A path to a directory of user own genome. The directory
#'      must contain chromosome-separated fasta files. The name of the fasta
#'      files must be similar to chromosome name in the coordinate table.
#' @field rm.case.kmer.overlaps Remove overlapping case k-mers? This is useful
#'      to completely obliterate potential neighboring effects.
#' @field single.case.len If the case is in uniform length, only start positions
#'      are considered. This is to remove redundancy of the coordinates and
#'      reduce memory usage.
#' @field merge.replicates Merge the replicates if exist? Merging will remove
#'      all duplicate positions from multiple replicates, resulting in a single
#'      unique individual point. Not merging replicate will result in weighted
#'      k-mer counts. Default is TRUE.
#' @field kmer.table K-mer table with pre-calculated scores. It can be directly
#'      used in other kmeRtone module.
#' @field module kmeRtone module. Options are: "score" (calculate z-score),
#'      "explore" (perform exploratory analysis), "tune" (find the best length
#'      of k-mer to be used).
#' @field rm.dup Remove duplicate coordinate? Remove duplicates within replicate
#'      e.g. sequencing read duplicates. Default is TRUE.
#' @field case.coor.1st.idx Case coordinate indexing format whether
#'      zero-based [start, end) or one-based [start, end]. Default is 1.
#' @field ctrl.coor.1st.idx Control coordinate indexing format whether
#'      zero-based [start, end) or one-based [start, end]. Default is 1.
#' @field coor.load.limit Maximum loaded chromosome coordinate data. Default is
#'      1.
#' @field genome.load.limit Maximum loaded chromosome sequence data. Default is
#'      1.
#' @field genome.fasta.style Genome FASTA style: "UCSC" or "NCBI". Default
#'      is "UCSC"
#' @field genome.ncbi.db For NCBI Genome, select database to use: "refseq" or
#'      "genbank". Default is "refseq".
#' @field use.UCSC.chr.name For NCBI Genome, use UCSC chromosome name?
#' @field verbose Print message? Default is TRUE.
#' @field kmer.cutoff Percent kmer cutoff for the case studies. Default is 5.
#' @field selected.extremophiles A vector of selected extremophile species. e.g.
#'    c("Deinococcus soli", "Deinococcus deserti")
#'    The best representative will be selected from the assembly summary.
#' @field other.extremophiles A vector of other extremophile species. These are
#'    used as a control to compare with the selected extremophiles.
#' @field cosmic.username COSMIC username (email). This is used to get cancer
#'    gene census.
#' @field cosmic.password COSMIC password.
#' @field tumour.type.regex Tumour-type regular expression to filter the cancer
#'    gene census table.
#' @field tumour.type.exact Exact tumour type to be included in the cancer gene
#'    census table.
#' @field cell.type Cell type to be included in the cancer gene census table.
#'    Default is somatic cell.
#' @field genic.elements.counts.dt Counts of susceptible k-mers in the genic
#'    elements, generated from STUDY_GENIC_ELEMENTS.
#' @field population.size Size of population in STUDY_ACROSS_POPULATION. Default
#'    is 1 million.
#' @field selected.genes Selected genes to be mutated based on SNVs in
#'    STUDY_ACROSS_POPULATION.
#' @field add.to.existing.population STUDY_ACROSS_POPULATION can be run multiple
#'    times to add to the existing simulated individuals. Default is FALSE to
#'    overwrite.
#' @field population.snv.dt SNVs used in the population simulation.
#' @field pop.plot To plot the outcome of STUDY_ACROSS_POPULATION. Default is
#'    TRUE.
#' @field pop.loop.chr In STUDY_ACROSS_POPULATION, should the operation loop
#'    based on chromosome name? Default is FALSE.

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
  cat("\n")

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
      cat("Argument single.case.len is not specified. Case length is assumed",
          " to be varied.\n")
      single.case.len <- NULL
    } else {
      single.case.len <- case$single_len
    }
  }
  if (missing(case.pattern) & module != "study_cancer_genes") {
    cat("Case pattern is not specified. Case pattern is assumed to be any",
        "pattern i.e. it is set to NULL.\n")
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