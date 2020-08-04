analyseGeneTranscripts <- function(genome.name, genome.annotation, k,
                                   ref.kmers, cutoff, ncpu,
                                   output = "data", igr.rel.pos = c(5000, 7500),
                                   igr.min.length = 150,
                                   gene.upstream.length = 2500,
                                   gene.downstream.length = 1000,
                                   tumour.keyword, tumour.exact.keyword,
                                   username = NULL, password = NULL,
                                   seed.number, replication.number,
                                   subsample.number = NULL,
                                   genome.path = NULL, genome.prefix = "",
                                   genome.suffix = NULL, genome = NULL,
                                   genome.annotation.path = NULL,
                                   transcript.susceptibility = NULL,
                                   cancer.gene.census = NULL
                                   ) {

  # AIM: To analyse human genome annotation kmeric patterns. The genome region
  #      is divided into exon, intron, 5' UTR, ...
  #
  # Workflow:
  # 1. Load genome annotation table.
  #      - Automatically download the table if it is not available locally.
  #        Default location data/genome/human/annotation
  #      - Existing table can

  # DNA pattern is inferred from the ref.kmers.

  # genome.name         <string>   Name of available human genome assembly:
  #                                "hg19" or "hg38".
  # genome.annotation   <string>   Genome annotation table name based on UCSC
  #                                Genome Browser. e.g. "ncbiRefSeq", "refGene",
  #                                "knownGene", etc.
  #                                ref: https://genome.ucsc.edu/cgi-bin/hgTables
  #                   <data.table> Pre-loaded table is acceptable but must
  #                                conform to the UCSC format.
  # k                   <numeric>  The size of a kmer.
  # ref.kmers           <string>   A path to scored kmers.
  #                   <data.table> A <data.table> format is accepted.

  # Second part of operation
  # Source of data (cancer_gene_census.csv) is downloded from
  # https://cancer.sanger.ac.uk/census
  # Abbreviation and other info can be found from that site as well.
  # The file is free for academic research communities.

  # transcript.           <string>    A path to transcript.susceptibility table
  #    susceptibility   <data.table>  or <data.table>. If NULL, the first
  #                                   operation will be skipped.
  # cancer.gene.census  <data.table>  A data.table of cancer gene census
  #                       <string>    A path to cancer gene census file
  # username              <string>    A username i.e. email registered to
  #                                   Sanger Cosmic
  # password              <string>    A password
  #
  # For boolean information, the table use yes and empty space for no.
  #
  # Dependencies
  #    Libraries: data.table, stringi

  env = environment()

  ## Dependant libraries #######################################################
  # order is important due to namespace masking effect

  libraries = c("stringi", "data.table", "base64enc", "R.utils")

  for (lib in libraries) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }

  ## Dependant functions #######################################################
  source("lib/loadGenomeV3.R")
  source("lib/splitFasta.R")
  source("lib/mergeGenCoordinate.R")
  source("lib/trimGenCoordinates.R")
  source("lib/removeCaseZone.R")
  source("lib/extractKmers.R")
  source("lib/reverseComplement.R")
  source("lib/distributeChunk.R")
  source("lib/distributeChunk2.R")
  source("lib/countReverseComplementKmers.R")
  source("lib/gene_transcript/processUCSCannoTable.R")
  source("lib/gene_transcript/calculateTranscriptSusceptibility.R")
  #source("lib/gene_transcript/drawGene.R")

  # Task specific dependant functions
  base.lib = "lib/gene_transcript/"
  source(paste0(base.lib, "01_checkInput.R"), local = TRUE)
  source("lib/02_prepGenome.R", local = TRUE)
  source(paste0(base.lib, "03_prepGenomeAnnotation.R"), local = TRUE)
  source(paste0(base.lib, "03b_downloadHumanGenomeAnnotation.R"), local = TRUE)
  source(paste0(base.lib, "04_processGenomeAnnotationTable.R"), local = TRUE)
  source(paste0(base.lib, "05_calculateSusceptibility.R"), local = TRUE)
  source(paste0(base.lib, "06_plotSusceptibility.R"), local = TRUE)
  source(paste0(base.lib, "07_downloadCancerGeneCensus.R"), local = TRUE)
  source(paste0(base.lib, "08_plotCancerGeneSusceptibilityDensity.R"))
  source(paste0(base.lib, "keyword2regex.R"))

  ## Parallel setup ############################################################
  if (ncpu > 1) {

    cpu.libs = c("foreach", "doParallel")
    for (lib in cpu.libs) {
      suppressPackageStartupMessages(library(lib, character.only = TRUE))
    }

    cl <- makeCluster(ncpu, outfile = "")
    registerDoParallel(cl)

  }

  ## Directory setup ###########################################################
  suppressWarnings(dir.create(output))

  if (is.null(transcript.susceptibility)) {

  # ---------------- A. INPUT CHECKING -----------------------------------------

  cat("[1] Checking inputs...")
  checkInput(genome.name, genome.annotation, k, ref.kmers, output.file,
             genome.path, genome.prefix, genome.suffix, genome.annotation.path)
  cat("DONE!\n")

  # ---------------- B. GENOME -------------------------------------------------
  # Load genome

  cat("[2] Loading genome...")
  prepGenome(genome.name, genome.path, genome.prefix, genome.suffix, genome,
             env)
  cat("DONE!\n")

  # ---------------- C. GENOME ANNOTATION --------------------------------------
  # Load genome annotation

  # Load genome annotation
  cat("[3] Loading genome annotation table...")
  prepGenomeAnnotation(genome.annotation, genome.name, genome.annotation.path,
                       env)
  cat("DONE!\n")

  # Extract transcript elements
  # Output: element.coordinate
  # If exist, skip
  cat("[4] Processing genome annotation table...")
  processGenomeAnnotationTable(genome.annotation, genome, igr.rel.pos,
                               igr.min.length, gene.upstream.length,
                               gene.downstream.length)
  fwrite(element.coordinate, paste0(output, "/",
                                    "transcript_element_coordinate.csv"))
  cat("DONE!\n")

  # ---------------- D.  -----------------------------------------
  # Count susceptibility percentage of every transcript elements
  # Output: transcript.susceptibility

  ref.kmers <- fread(ref.kmers)

  cat("[5] Calculating transcript susceptibility...\n")
  calculateSusceptibility(element.coordinate, genome, ref.kmers, cutoff,
                          output)
  cat("DONE!\n")

  # Skip the step above if user input transcript.susceptibility table
  } else if (!is.null(transcript.susceptibility)) {
    cat("Skipping kmer susceptibility counting...\n")

    if (is.character(transcript.susceptibility)) {
      transcript.susceptibility <- fread(transcript.susceptibility)
    } else if (!is.data.table(transcript.susceptibility)) {
      stop("Please input transcript.susceptibility as file path or data.table")
    }
  }

  cat("[6] Plotting boxplot...")
  #plotSusceptibility(transcript.susceptibility, output)
  cat("DONE!\n")

  # ---------------- E. Prepare Cancer Gene Census -----------------------------
  # Download cancer_gene_census.csv from COSMIC sanger website

  if (is.character(cancer.gene.census)) {
    cancer.gene.census <- fread(cancer.gene.census, na.strings="")
  } else if (!(is.null(username) & is.null(password))) {
      cat("[7] Downloading cancer_gene_census.csv...")
      downloadCancerGeneCensus(username, password, output)
      cancer.gene.census <- paste0(output, "/cancer_gene_census.csv")
      cat("DONE!\n")
  }



  # ---------------- F. Cancer Gene Susceptibility -----------------------------
  # Calculate and plot density of cancer gene susceptibility

  cat("[8] Calculating density of cancer gene susceptibility...")
 # plotCancerGeneSusceptibilityDensity(cancer.gene.census,
 #                                     transcript.susceptibility,
 #                                     tumour.exact.keyword,
 #                                     tumour.regex.keyword,
 #                                     output)
  source("../kmertone/lib/gene_transcript/bootstrappingGenes.R")
bootstrappingGenes(transcript.susceptibility, cancer.gene.census,
                              tumour.keyword, tumour.exact.keyword,
                              seed.number, replication.number,
                              subsample.number, output, ncpu)
  cat("DONE!\n")
  # ---------------- THE END ---------------------------------------------------

  if (ncpu > 1) {
    stopCluster(cl)
  }

  return(transcript.susceptibility) # kmers table: kmer, count, fold_change, p, z
}
