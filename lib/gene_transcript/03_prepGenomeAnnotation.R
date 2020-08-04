prepGenomeAnnotation <- function(genome.annotation, genome.name,
                                 genome.annotation.path=NULL,
                                 env=parent.frame()) {
  # AIM: To load genome annotation table
  # Objectives: (1) To check whether the table exist locally
  #                 genome.annotation.path
  #             (2) To download the table if not exist.
  #             (3) To load the table.

  # If the user input data.table of genome annotation table
  if (class(genome.annotation)[1] == "data.table")
    return(cat("\n -- The table is already loaded.\n"))

  if (class(genome.annotation)[1] != "character")
    stop("Please input genome.annotation as <string> or <data.table>.")

  if (is.null(genome.annotation.path)) {
    genome.annotation.path <- paste0("data/genome/human/", genome.name, "/",
                                     genome.annotation, ".txt.gz")
  }

  if (!file.exists(genome.annotation.path)) {
    cat("\n")
    downloadHumanGenomeAnnotation(genome.name, genome.annotation,
                                  output.dir=paste0("data/genome/human/",
                                                    genome.name))
  }

  # For NCBI RefSeq or UCSC RefSeq annotation
  # The table has 16 columns. For further clarification, please refer to
  # http://genome.ucsc.edu/cgi-bin/hgTables, specifically click "describe tabl
  # schema" next to "table" selection.

  if (genome.annotation %in% c("ncbiRefSeq", "refGene")) {
    gen.annot.columns <- c("bin", "name", "chrom", "strand", "txStart", "txEnd",
                           "cdsStart",	"cdsEnd", "exonCount", "exonStarts",
                           "exonEnds", "score",	"name2", "cdsStartStat",
                           "cdsEndStat", "exonFrames")

    # For GENCODE annotation
  } else if (genome.annotation %in% c("knownGene")) {
    gen.annot.columns <- c("name", "chrom", "strand", "txStart", "txEnd",
                           "cdsStart", "cdsEnd", "exonCount", "exonStarts",
                           "exonEnds", "proteinID", "alignID")
  } else {
    stop(paste0("Sorry.", genome.annotation, "is not supported. Please message",
                "developer to support it. Thank you."))
  }

  assign("genome.annotation", fread(genome.annotation.path,
                                    col.names = gen.annot.columns,
                                    showProgress = FALSE),
         envir = env)

}
