analyseHumanTranscripts <- function(genome.name, genome.annotation, k, ref.kmers, score, cutoff, ncpu,
                                    output="data", igr.rel.pos = c(5000, 7500), igr.min.length = 150,
                                    gene.upstream.length = 2500, gene.downstream.length = 1000, 
                         genome.path=NULL, genome.prefix="", genome.suffix=NULL, genome=NULL,
                         genome.annotation.path=NULL) {
  
  # AIM: To analyse human genome annotation kmeric patterns. The genome region is divided into exon, intron,
  #      5' UTR, ...
  # Workflow:
  # 1. Load genome annotation table.
  #      - Automatically download the table if it is not available locally. Default location data/genome/human/annotation
  #      - Existing table can 
  
  # DNA pattern is inferred from the ref.kmers.
  
  # genome.name                  <string>         Name of available human genome assembly: "hg19" or "hg38").
  # genome.annotation            <string>         Genome annotation table name based on UCSC Genome Browser.
  #                                               e.g. "ncbiRefSeq", "refGene", "knownGene", etc.
  #                                               ref: https://genome.ucsc.edu/cgi-bin/hgTables
  #                             <data.table>      Pre-loaded table is acceptable but must conform to the UCSC format.
  # k                            <numeric>        The size of a kmer.
  # ref.kmers                    <string>         A path to scored kmers.
  #                             <data.table>      A <data.table> format is accepted.

  igr.rel.pos = c(5000, 7500)
  igr.min.length = 150
  gene.upstream.length = 2500
  gene.downstream.length = 1000
  score = "z"
  cutoff=95
  ref.kmers <- "data/kmers_k-8_pattern-CT.csv"
  output="data"
  ncpu=6
  
  genome.name="hg38"
  genome.annotation="ncbiRefSeq"
  output.folder="output"
  #DNA.pattern=c("TT")
  k=8
  
  genome.path=NULL; genome.prefix=""; genome.suffix=NULL; genome=NULL; genome.annotation.path=NULL
  
  env = environment()
  #setwd("../kmertone/")
  
  # location of the TrantoR library
  #TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant libraries #########################################################
  # order is important due to namespace masking effect
  #suppressPackageStartupMessages( library(Biostrings) )
  suppressPackageStartupMessages( library(  stringi ) )
  suppressPackageStartupMessages( library(data.table) )
  
  ## Dependant functions #########################################################
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
  source("lib/Human_genome_annotation/processUCSCannoTable.R")
  source("lib/Human_genome_annotation/calculateTranscriptSusceptibility.R")
  #source("lib/Human_genome_annotation/drawGene.R")

  # Dependant functions from the TrantorR library
  #source("lib/TrantoRext/GEN_getSeqLogo.R", local = TRUE) # dependent on seqlogo from bioconductor
  
  # Task specific dependant functions
  base.lib = "lib/Human_genome_annotation/"
  source(paste0(base.lib, "01_checkInput.R"), local = TRUE)
  source("lib/02_prepGenome.R", local = TRUE)
  source(paste0(base.lib, "03_prepGenomeAnnotation.R"), local = TRUE)
  source(paste0(base.lib, "03b_downloadHumanGenomeAnnotation.R"), local = TRUE)
  source(paste0(base.lib, "04_processGenomeAnnotationTable.R"), local = TRUE)
  source(paste0(base.lib, "05_calculateSusceptibility.R"), local = TRUE)
  
  #source("lib/07_getCaseKmers.R", local = TRUE)
  #source("lib/08_getControlKmers.R", local = TRUE)
  #source("lib/XX_mapZscore.R", local = TRUE)

  ## Parallel setup ##############################################################
  if (ncpu > 1) {

    suppressPackageStartupMessages( library(foreach)    )
    suppressPackageStartupMessages( library(doParallel) )

    cl <- makeCluster(ncpu, outfile="")
    registerDoParallel(cl)

  }
  
  ## Directory setup #############################################################
  suppressWarnings(dir.create(output))
  
  
  # ---------------- A. INPUT CHECKING -----------------------------------------------------
  
  cat("[1] Checking inputs...")
  checkInput(genome.name, genome.annotation, k, ref.kmers, output.file,
             genome.path, genome.prefix, genome.suffix, genome.annotation.path)
  cat("DONE!\n")
  
  # ---------------- B. GENOME -------------------------------------------------------------
  # Load genome
  
  cat("[2] Loading genome...")
  prepGenome(genome.name, genome.path, genome.prefix, genome.suffix, genome, env)
  cat("DONE!\n")
  
  # ---------------- C. GENOME ANNOTATION --------------------------------------------------
  # Load genome annotation
  
  # Load genome annotation
  cat("[3] Loading genome annotation table...")
  prepGenomeAnnotation(genome.annotation, genome.name, genome.annotation.path, env)
  cat("DONE!\n")
  
  # Extract transcript elements
  # Output: element.coordinate
  cat("[4] Processing genome annotation table...")
  processGenomeAnnotationTable(genome.annotation, genome, igr.rel.pos, igr.min.length,
                               gene.upstream.length, gene.downstream.length)
  fwrite(element.coordinate, paste0(output.folder, "/", "transcript_element_coordinate.csv"))
  cat("DONE!\n")
  
  # ---------------- D. SUSCEPTIBILITY -----------------------------------------------------
  # Count susceptibility percentage of every transcript elements
  # Output: transcript.susceptibility
  
  cat("[5] Calculating transcript susceptibility...")
  calculateSusceptibility(element.coordinate, genome, ref.kmers, score, cutoff, output)
  cat("DONE!\n")
  
  cat("[6] Plotting boxplot...")
  plotSusceptibility(susceptibility.score, ouput)
  cat("DONE!\n")
  
  # ---------------- E. KMER SCORE MAPPING ------------------------------------------------
  # Map Z score to the kmers
  
  cat("[5] Mapping kmer scores...")
  mapZscores()
  cat("DONE!\n")
  
  # ---------------- PRE-ANALYSIS --------------------------------------------------------
  # Replicate summary
  # GC and G content at various width
  
  # case distribution
  #cat("Distribution of case site after filtering.\n\n")
  #caseDistribution("genomic.coordinate", DNA.pattern, output.file=output.file, env=kmertone.env)
  
  # seqlogos
  
  
  # Optimum K
  #calculateOptimumK(kmertone.env, set.k = c(2,4,6,8,10,12))
  
  # G|C and G content
  #GCcontent(env)
  
  # Plot G|C and G density
  #plotDensity() # TBD
  
  # seqlogo
  #drawSeqlogo() # TBD
  
  # volcano plot
  #plotVolcano() # TBD
  
  # ---------------- KMER EXTRACTION ------------------------------------------------------
  
  # cat("[5] Getting case kmers...\n\n")
  # kmers <- getCaseKmers("genomic.coordinate", genome, k, DNA.pattern, strand.mode,
  #              remove.overlaps=TRUE, kmertone.env)
  # fwrite(kmers, paste0(output.file, "/kmers.csv"))
  # cat("\n")
  # 
  # if (!is.null(control.relative.position)) {
  #   cat("[6] Getting control kmers...\n\n")
  #   kmers <- getControlKmers("genomic.coordinate", genome, k, DNA.pattern, strand.mode,
  #                            "kmers", kmertone.env)
  #   fwrite(kmers, paste0(output.file, "/kmers.csv"))
  # }
  # 
  # # ---------------- UPDATE PRE-ANALYSIS ---------------------------------------------------
  # 
  # #GCcontent(dts[1], genome, filename = "GC_after")
  # #Gcontent(dts[1], genome, filename = "G_after")
  # #seqlogo(dts[1], genome, filename = "seqlogo_after")
  # 
  # # ---------------- SCORE -----------------------------------------------------------------
  # 
  # cat("\n[7] Calculating z score...")
  # zScore("kmers", kmertone.env)
  # #pValue("kmers") # TBD
  # cat("DONE!\n")
  
  # ---------------- THE END ---------------------------------------------------------------
  
  if (ncpu > 1) {
    stopCluster(cl)
  }
  
  #return(kmers) # kmers table: kmer, count, fold_change, p, z
}