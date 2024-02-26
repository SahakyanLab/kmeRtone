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
    control <- buildControl(case = case,
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
