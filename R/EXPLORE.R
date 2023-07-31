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
