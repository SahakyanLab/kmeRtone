#' Calculate susceptibility scores for k-mers in case and control regions.
#'
#' Function calculates susceptibility scores for k-mers in case and control regions.
#' Case regions are defined by genomic coordinates provided in a file or data.table.
#' Control regions can be constructed relative to the case regions or provided directly.
#' The scores are computed based on the occurrence of k-mers in case and control regions.
#'
#' @param case.coor.path Path to the file containing genomic coordinates of case regions.
#' @param genome.name Name of the genome to be used.
#' @param strand.sensitive Logical indicating whether strand information should be considered.
#' @param k Length of the k-mers to be extracted.
#' @param ctrl.rel.pos Relative positions of control regions with respect to case regions. 
#'                      It should be a vector of two integers indicating the upstream and downstream 
#'                      distances from the case regions.
#' @param case.pattern Regular expression pattern to identify the central sequence in case regions.
#' @param output.path Directory path where the output files will be saved.
#' @param case Data.table containing the genomic coordinates of case regions.
#' @param genome Genome data.table containing the genomic sequence information.
#' @param control Data.table containing the genomic coordinates of control regions.
#' @param control.path Path to the file containing genomic coordinates of control regions (optional).
#' @param genome.path Path to the genome FASTA file.
#' @param rm.case.kmer.overlaps Logical indicating whether overlapping k-mers within case regions should be removed.
#' @param single.case.len Single case length.
#' @param merge.replicates Logical indicating whether replicates should be merged.
#' @param rm.dup Logical indicating whether duplicate k-mers should be removed.
#' @param case.coor.1st.idx First index in the case coordinate file.
#' @param ctrl.coor.1st.idx First index in the control coordinate file.
#' @param coor.load.limit Maximum number of coordinates to load.
#' @param genome.load.limit Maximum number of genome sequences to load.
#' @param genome.fasta.style FASTA style.
#' @param genome.ncbi.db NCBI database.
#' @param use.UCSC.chr.name Logical indicating whether to use UCSC chromosome names.
#' @param verbose Logical indicating whether to display progress messages.
#'
#' @importFrom data.table fread fwrite
#' @importFrom stringi stri_replace_first_fixed stri_replace_all_regex
#' stri_split_fixed stri_sub stri_locate_first_regex stri_replace_first_regex
#' stri_match_all_regex
#'
#' @return Data.table containing susceptibility scores for k-mers.
#'
#' @export
SCORE <- function(
  case.coor.path, genome.name, strand.sensitive, k, ctrl.rel.pos, case.pattern,
  output.path, case, genome, control, control.path, genome.path,
  rm.case.kmer.overlaps, single.case.len, merge.replicates,
  rm.dup, case.coor.1st.idx, ctrl.coor.1st.idx, coor.load.limit,
  genome.load.limit, genome.fasta.style, genome.ncbi.db, use.UCSC.chr.name,
  verbose) {

  if (is.null(case))
    case <- loadCoordinate(root.path = case.coor.path,
                           single.len = single.case.len,
                           is.strand.sensitive = strand.sensitive,
                           merge.replicates = merge.replicates,
                           rm.dup = rm.dup,
                           add.col.rep = FALSE,
                           is.kmer = FALSE,
                           ori.first.index = case.coor.1st.idx,
                           load.limit = coor.load.limit)

  if (is.null(genome))
    genome <- loadGenome(fasta.path = genome.path,
                         genome.name = genome.name,
                         fasta.style = genome.fasta.style,
                         ncbi.db = genome.ncbi.db,
                         mask = "none",
                         use.UCSC.name = use.UCSC.chr.name,
                         load.limit = genome.load.limit)

  # Time
  if (verbose) T1 <- t1 <- Sys.time()

  if (verbose) catHeader("Extraction of Case K-mers")
  case.kmers <- extractKmers(coor = case,
                             genome = genome,
                             k = k,
                             central.pattern = case.pattern,
                             rm.overlap.region = rm.case.kmer.overlaps,
                             verbose = verbose)

  # Time
  if (verbose) {
    t <- Sys.time() - t1
    cat("\nTotal time taken:", round(t[[1]], 2), attr(t, "units"),
        "\n")
    t1 <- Sys.time()
  }

  if (verbose) catHeader("Extraction of Control K-mers")
  if (is.null(control) & is.null(control.path)) {
  #' @param case Data.table containing the genomic coordinates of case regions.
  #' @param k Integer size of the expanded k-mer.
  #' @param ctrl.rel.pos Relative positions of control regions with respect to case regions. 
  #'    It should be a vector of two integers indicating the upstream and downstream
  #'    distances from the case regions.
  #' @param genome Genome data.table containing the genomic sequence information.
  #' @param output.path Directory path where the output files will be saved.
  #' @param verbose Logical indicating whether to display progress messages.
    control <- buildControl(case = case,
                            k = k,
                            ctrl.rel.pos = ctrl.rel.pos,
                            genome = genome,
                            output.path = paste0(output.path, "/control_",
                                                 ctrl.rel.pos[1], "-",
                                                 ctrl.rel.pos[2], "/"),
                            verbose = verbose)

    # Time
    if (verbose) {
      t <- Sys.time() - t1
      cat("\nTotal time taken:", round(t[[1]], 2), attr(t, "units"), "\n")
      t1 <- Sys.time()
    }
  } else if (!is.null(control.path)) {
    control <- Coordinate$new(root.path = control.path,
                              is.strand.sensitive = FALSE,
                              ori.first.index = ctrl.coor.1st.idx)
  }

  control.kmers <- extractKmers(coor = control,
                                genome = genome,
                                k = k,
                                central.pattern = case.pattern,
                                rm.overlap.region = FALSE,
                                verbose = verbose)

  # Time
  if (verbose) {
    t <- Sys.time() - t1
    cat("\nTotal time taken:", round(t[[1]], 2), attr(t, "units"), "\n")
    t1 <- Sys.time()
  }

  if (verbose) catHeader("Calculation of K-mer Susceptibility")
  kmer.table <- getScores(case.kmers = case.kmers,
                          control.kmers = control.kmers)

  fwrite(kmer.table, paste0(output.path, "/score_", k, "-mers.csv"),
         showProgress = FALSE)

  message("The ", k, "-mer scores are saved at ", output.path, "/score_",
          k,"-mer.csv")

  # Time
  if (verbose) {
    t <- Sys.time() - T1
    cat("\nFINISH! Total time taken:", round(t[[1]], 2), attr(t, "units"),
        "\n")
  }

  return(kmer.table)
}
