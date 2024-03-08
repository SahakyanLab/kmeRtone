#' Resolve and generate genic element coordinates from UCSC genePred table.
#'
#' Function generates intergenic coordinates from a UCSC genePred table. 
#' It allows users to specify the genePred data source, the relative position 
#' and minimum length for intergenic regions, and whether to return the results 
#' as a `Coordinate` object or a `data.table`.
#'
#' @param genepred UCSC genePred table or database name ("refseq" or "gencode").
#' @param genome.name UCSC genome name (e.g., hg38, mm39).
#' @param igr.rel.pos Intergenic relative position, defaults to c(5000, 7500).
#' @param igr.min.length Minimum length for intergenic regions, default is 150.
#' @param return.coor.obj Return results as a `Coordinate` object? Default FALSE.
#' @return Intergenic coordinates as a `data.table` or `Coordinate` object.
#'
#' @importFrom data.table data.table fwrite
#' 
#' @export
generateIntergenicCoor <- function(genepred, genome.name,
                                   igr.rel.pos=c(5000, 7500),
                                   igr.min.length=150, return.coor.obj=FALSE) {

  genome <- loadGenome(genome.name, fasta.style = "UCSC")

  if (is.character(genepred))
    genepred <- getUCSCgenePredTable(genome.name = genome.name, db = genepred)

  # Filter out non-protein-coding genes.
  col.name <- c("transcriptClass", "name")
  col.name <- col.name[col.name %in% names(genepred)][1]
  genepred <- genepred[get(col.name) %like% "^NM_|^XM_|coding"]

  # Create gene regions.
  gene <- genepred[, .(chromosome = chrom, start = txStart + 1, end = txEnd)] |>
    mergeCoordinate()

  # Create intergenic regions.
  igr <- gene[, .(chromosome,
                  start = c(start - igr.rel.pos[2],
                            end + igr.rel.pos[1]),
                  end = c(start - igr.rel.pos[1],
                          end + igr.rel.pos[2]))] |>
    trimCoordinate(genome = genome) |>
    mergeCoordinate()

  ## Add buffer padding to the gene regions.
  gene[, `:=`(start = start - igr.rel.pos[1],
              end = end + igr.rel.pos[1])]

  # Remove gene region + buffer from the intergenic region.
  igr <- removeRegion(coor = igr, region = gene)

  # Restrict intergenic size.
  igr <- igr[end - start + 1 >= igr.min.length]

  # Name the intergenic regions.
  igr[, element := "IGR"]

  if (return.coor.obj) {
    tmp.dir <- tempfile()
    igr[, fwrite(.SD, paste0(tmp.dir, "/", chromosome, ".csv"),
                 showProgress = FALSE),
        by = chromosome]
    igr <- loadCoordinate(root.path = tmp.dir,
                          single.len = NULL,
                          is.strand.sensitive = TRUE,
                          merge.replicates = FALSE,
                          rm.dup = FALSE,
                          add.col.rep = FALSE,
                          is.kmer = FALSE,
                          k = NA,
                          ori.first.index = 1,
                          load.limit = 1)
  }

  return(igr)
}
