kmertone <- function(genomic.coordinate, genome.name, strand.sensitive, k,
                     control.relative.position, DNA.pattern=NULL, output="data",
                     genome=NULL, genome.path=NULL, genome.prefix="",
                     genome.suffix=NULL, shrink.coordinate=FALSE, ncpu=1) {

  # Kmertone program follows UCSC Genome convention e.g. chromosomes are named
  # as chr1, chr2, chr3, chr?, ...
  # However indexing format strictly follows R format i.e. one-based index.
  # Genomic coordinate table is a BED-like format with only 4 columns in order:
  # chromosome, start, end, strand.

  # Flag                   Format     Description
  # genomic.coordinate  <data.table>  A 0-based index genomic coordinate. A list
  #                                   of data.table denotes multiple replicates.
  #                                   The table columns must be (in order):
  #                                   chromosome, start, end, and strand.
  #                       <string>    A path to genomic coordinate. A vector of
  #                                   paths denotes multiple replicates. BED
  #                                   file is supported.
  #
  # genome.name           <string>    Name of available genome: "hg19" or "hg38"
  #
  # genome.path           <string>    A folder path to user own genome. The
  #                                   folder must contain separated chromosome
  #                                   fasta files. The name of the fasta files
  #                                   should be similar to chromosome name in
  #                                   the genomic coordinate table.
  #
  # genome.suffix         <string>    An extension name of the fasta files.
  #                                   Compressed file is supported.
  #
  # DNA.pattern           <string>    A single or multiple DNA pattern.
  #                                   e.g. "G", "TT", "TC", etc.
  #                                   It can be set as NULL.
  #
  # k                     <numeric>   The size of a kmer.
  #                       <string>    User can input "optimum" to calculate and
  #                                   use optimum k.
  #
  # control.relative.     <vector>    A relative position of control regions in
  #    position                       a vector format i.e. c(start, end).
  #                                   For example c(80,500) means control
  #                                   regions are 80-500 bases away from the
  #                                   the case site (upstream and downstream)
  #
  # strand.sensitive       <bool>     Mode of strand. In sensitive mode, case
  #                                   kmers are extracted from the sense strand
  #                                   only. In contrast, in insensitive mode,
  #                                   the case kmers are extracted from both.
  #
  # shrink.coordinate      <bool>     To shrink genomic coordinate down to case
  #                                   coordinate only.
  #
  # ncpu                  <numeric>   Number of cpu core to use. Only one core
  #                                   is supported at the moment.
  ncpu=1

  kmertone.env = environment()
  #setwd("../kmertone/")

  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"

  ## Dependant libraries #######################################################
  # order is important due to namespace masking effect
  suppressPackageStartupMessages( library(Biostrings)  )
  suppressPackageStartupMessages( library(  seqLogo )  )
  suppressPackageStartupMessages( library( venneuler)  )
  suppressPackageStartupMessages( library(  stringi )  )
  suppressPackageStartupMessages( library(data.table)  )
  suppressPackageStartupMessages( library(RColorBrewer))

  ## Dependant functions #######################################################
  source("lib/loadGenomeV3.R")
  source("lib/splitFasta.R")
  source("lib/reverseComplement.R")
  source("lib/distributeChunk.R")
  source("lib/distributeChunk2.R")
  source("lib/scaleGenCoordinate.R")
  source("lib/trimGenCoordinates.R")
  source("lib/mergeGenCoordinate.R")
  source("lib/expandGenCoordinate.R")
  source("lib/extractKmers.R")
  source("lib/countReverseComplementKmers.R")
  source("lib/removeAllOverlaps.R")
  source("lib/removeCaseZone.R")
  source("lib/extractGenomeKmers.R")

  # Dependant functions from the TrantorR library
  source(paste0(TrantoRLib, "GEN_getSeqLogo.R"), local = TRUE)

  # Task specific dependant functions
  source("lib/01_inputChecking.R", local = TRUE)
  source("lib/02_prepGenome.R", local = TRUE)
  source("lib/03a_prepGenCoordinate.R", local = TRUE)
  source("lib/03b_addColumnSequence.R", local = TRUE)
  source("lib/03c_filterTable.R", local = TRUE)
  source("lib/04_caseDistribution.R", local = TRUE)
  #source("lib/05_GCcontent.R", local = TRUE)
  source("lib/07_getCaseKmers.R")
  source("lib/08_getControlKmers.R", local = TRUE)
  source("lib/08_zScore.R", local = TRUE)
  source("lib/05_calculateOptimumK.R")
  #source("lib/seqlogo.R", local = TRUE)

  ## Parallel setup ############################################################
  if (ncpu > 1) {
    suppressPackageStartupMessages( library(foreach)    )
    suppressPackageStartupMessages( library(doParallel) )
    cl <- makeCluster(ncpu)
    registerDoParallel(cl)
  }

  ## Directory setup ###########################################################
  suppressWarnings(dir.create(output, recursive = TRUE))


  # ---------------- A. INPUT CHECKING -----------------------------------------
  cat("[1] Checking inputs...")
  inputChecking(DNA.pattern, k, control.relative.position, strand.sensitive)
  cat("DONE!\n")

  # ---------------- B. GENOME -------------------------------------------------
  # 1. Load genome

  cat("[2] Loading genome...")
  prepGenome(genome.name, genome.path, genome.prefix, genome.suffix, genome,
             kmertone.env)
  cat("\n")

  # ---------------- C. GENOMIC COORDINATE -------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)

  cat("[3] Loading genomic coordinate table...\n")
  prepGenCoordinate("genomic.coordinate", strand.sensitive, genome,
		    kmertone.env)
  gc()
  cat("\n")

  # add column sequence if strand sensitive
  if (strand.sensitive){
    cat("[3a] Filtering genomic coordinate table...\n")
    addColumnSequence("genomic.coordinate", genome, shrink.coordinate,
                      kmertone.env)
    filterTable("genomic.coordinate", DNA.pattern, strand.sensitive, output,
                kmertone.env) # memory spike here
    fwrite(genomic.coordinate, paste0(output, "/filtered_table.csv"))
  }

  # Backup original coordinates
  genomic.coordinate[, c("original_start", "original_end") := list(start, end)]

  # ---------------- 4. PRE-ANALYSIS -------------------------------------------
  # Replicate summary
  # GC and G content at various width

  # case distribution
  cat("Distribution of case site after filtering.\n\n")
  caseDistribution("genomic.coordinate", DNA.pattern, output, env=kmertone.env)
  cat("\n")

  # Optimum K
  if(k == "optimum"){
    cat("Calculating optimum k...\n")

    q <- calculateOptimumK(genomic.coordinate, genome, set.k, DNA.pattern,
                      strand.sensitive, env=kmertone.env)

    cat("Optimum k is", k)
    # For testing purpose
    return(q)
  }


  # seqlogos

  # G|C and G content
  #GCcontent(env)

  # Plot G|C and G density
  #plotDensity() # TBD

  # seqlogo
  #drawSeqlogo() # TBD

  # volcano plot
  #plotVolcano() # TBD

  # ---------------- 5. KMER EXTRACTION ----------------------------------------

  cat("[5] Getting case kmers...\n\n")
  kmers <- getCaseKmers("genomic.coordinate", genome, k, DNA.pattern,
			strand.sensitive, remove.overlaps=TRUE, kmertone.env)
  cat("\n")

  genomic.coordinate[, `:=`(start = original_start, end = original_end)]
  if (!is.null(control.relative.position)) {
    cat("[6] Getting control kmers...\n\n")
    kmers <- getControlKmers("genomic.coordinate", genome,
			     control.relative.position, k, DNA.pattern,
			     strand.sensitive, "kmers", kmertone.env)
  }

  # # ---------------- 6. UPDATE PRE-ANALYSIS ----------------------------------
  #
  # #GCcontent(dts[1], genome, filename = "GC_after")
  # #Gcontent(dts[1], genome, filename = "G_after")
  # #seqlogo(dts[1], genome, filename = "seqlogo_after")
  #
  # # ---------------- 7. SCORE ------------------------------------------------

  cat("\n[7] Calculating z score...")
  zScore("kmers", DNA.pattern, kmertone.env)
  #pValue("kmers") # TBD
  cat("DONE!\n")

  # ---------------- THE END ---------------------------------------------------

  # if (ncpu > 1) {
  #   stopCluster(cl)
  # }

  # Save kmers
  if(is.null(DNA.pattern)){
    fwrite(kmers, paste0(output, "/kmers_k-", k, ".csv"))
  } else {
    fwrite(kmers, paste0(output, "/kmers_k-", k, "_pattern-",
			 paste(DNA.pattern, collapse = "_"), ".csv"))
  }

  return(kmers) # kmers table: kmer, count, fold_change, p, z
}
