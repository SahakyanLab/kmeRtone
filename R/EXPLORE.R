#' Function generates various exploratory analyses.
#'
#' @param case.coor.path Path to case coordinates.
#' @param genome.name Genome name (e.g., hg19, hg38).
#' @param strand.sensitive Boolean indicating if strand sensitivity is considered.
#' @param k K-mer size.
#' @param case.pattern String patterns to consider in the analysis.
#' @param output.path Output directory path for exploration plots.
#' @param case Coordinate class object or similar structure for case data.
#' @param genome Genome class object or similar structure.
#' @param control Control class object or similar structure.
#' @param genome.path Path to genome fasta files.
#' @param single.case.len Length of single cases.
#' @param rm.dup Boolean indicating if duplicates should be removed.
#' @param case.coor.1st.idx Indexing of case coordinates.
#' @param coor.load.limit Maximum number of coordinates to load.
#' @param genome.load.limit Maximum number of genome data to load.
#' @param genome.fasta.style Fasta file style for genome data.
#' @param genome.ncbi.db NCBI database for genome data.
#' @param use.UCSC.chr.name Boolean indicating if UCSC chromosome naming is used.
#' @param verbose Boolean indicating if verbose output is enabled.
#'
#' @return Output directory containing exploration plots.
#'
#' @importFrom R6 R6Class
#' @importFrom grDevices cairo_pdf dev.off
#' 
#' @export
EXPLORE <- function(
  case.coor.path, genome.name, strand.sensitive, k, case.pattern, output.path,
  case, genome, control, genome.path, single.case.len, rm.dup,
  case.coor.1st.idx, coor.load.limit, genome.load.limit, genome.fasta.style,
  genome.ncbi.db, use.UCSC.chr.name, verbose) {

  dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(paste0(output.path, "/exploration.pdf"), width = 10, height = 8,
            onefile = TRUE)

  if (is.null(case))
    case <- loadCoordinate(root.path = case.coor.path,
                           single.len = single.case.len,
                           is.strand.sensitive = strand.sensitive,
                           merge.replicates = TRUE,
                           rm.dup = rm.dup,
                           add.col.rep = TRUE,
                           is.kmer = FALSE,
                           ori.first.index = case.coor.1st.idx,
                           load.limit = coor.load.limit)
  else case$add_col_rep <- TRUE

  if (is.null(genome))
    genome <- loadGenome(fasta.path = genome.path,
                         genome.name = genome.name,
                         fasta.style = genome.fasta.style,
                         ncbi.db = genome.ncbi.db,
                         mask = "none",
                         use.UCSC.name = use.UCSC.chr.name,
                         load.limit = genome.load.limit)

  countDistribution(case, genome, case.pattern, output.path = NULL)
  countBaseComposition(case, genome, case.pattern, output.path = NULL)

  dev.off()

  if (verbose) cat("The exploration plots are saved at",
                   paste0(output.path, "/exploration.pdf\n"))
}
