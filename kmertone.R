kmertone <- function(genomic.coordinate, genome.name, strand.mode, k,
                     control.relative.position, DNA.pattern=NULL, output="data",
                     genome, genome.path=NULL, genome.prefix="", genome.suffix=NULL) {
  
  # Kmertone program follows UCSC Genome convention e.g. chromosomes are named as chr1, chr2, chr3, chr?, ...
  # However indexing format strictly follows R format i.e. one-based index.
  # Genomic coordinate table is a BED-like format with only 4 columns in order: chromosome, start, end, strand.
  
  # Flag                          Format          Description
  # genomic.coordinate           <data.table>     A data.table of genomic coordinate. It can be a list of
  #                                               data.table for replicates. It must contains column 
  #                                               chromosome, start, end, and strand. For "insensitive" mode
  #                                               strand information is not required.
  # genome.name                  <string>         Name of available genome: "GRCh37" or "GRCh38"
  # genome.path                  <string>         If provided genome is not available, user can input a path
  #                                               to a folder containing chromosome fasta files. The fasta
  #                                               files should be named as chromosome names just like in the
  #                                               genomic coordinate table and has a uniform extension name.
  #                                               The extension can  be specified in flag suffix below. The
  #                                               fasta files are expected to contain only one-line header
  #                                               followed by same-length sequence in every line.
  #                                               e.g. chr1.fasta, chr2.fasta, etc. The 
  # genome.suffix                <string>         An extension name of the fasta files. Compressed file is
  #                                               supported. e.g. .fasta, .fa, .fa.gz, etc.
  # DNA.pattern                  <string>         A single or multiple DNA pattern. e.g. "G", "TT", "TC", etc.
  #                                               It can be set as NULL if no pattern is desired.
  # k                            <numeric>        The size of a kmer.
  #                              <string>         "optimum" - Calculate optimum k and use it.
  # control.relative.position    <vector>         A coordinate in a vector format i.e. c(start, end) pointing
  #                                               to where the control kmers should be extracted from. The
  #                                               coordinate is relative to the damage site (upstream and
  #                                               downstream)
  # strand.mode                  <string>         Mode of strand: "sensitive" or "insensitive".
  #                                               In "sensitive" mode, strand information is important. The
  #                                               kmers are extracted from the sense strand only.
  #                                               In contrast, in "insensitive" mode, the kmers are extracted
  #                                               from both plus and minus strands.

  # library(data.table)
  # genomic.coordinate <- c("data/table.csv")
  # genome.name="hg19";genome.path=NULL;genome.prefix=""; genome.suffix=".fa.gz"
  # DNA.pattern="TT"; k=10; control.relative.position=c(80,500)
  # strand.mode="sensitive"; ncpu=1

  kmertone.env = environment()
  #setwd("../kmertone/")
  
  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant libraries #########################################################
  # order is important due to namespace masking effect
  suppressPackageStartupMessages( library(Biostrings) )
  suppressPackageStartupMessages( library(  seqLogo ) )
  suppressPackageStartupMessages( library( venneuler) )
  suppressPackageStartupMessages( library(  stringi ) )
  suppressPackageStartupMessages( library(data.table) )
  
  ## Dependant functions #########################################################
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
  source("lib/TrantoRext/GEN_getSeqLogo.R", local = TRUE) # dependent on seqlogo from bioconductor
  
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

  ## Parallel setup ##############################################################
  # if (ncpu > 1) {
  #   
  #   suppressPackageStartupMessages( library(foreach)    )
  #   suppressPackageStartupMessages( library(doParallel) )
  #   
  #   cl <- makeCluster(ncpu)
  #   registerDoParallel(cl)
  # 
  # }
  
  ## Directory setup #############################################################
  suppressWarnings(dir.create(output, recursive = TRUE))

  
  # ---------------- A. INPUT CHECKING -----------------------------------------------------
  
  cat("[1] Checking inputs...")
  inputChecking(DNA.pattern, k, control.relative.position, strand.mode)
  cat("DONE!\n")
  
  # ---------------- B. GENOME -------------------------------------------------------------
  # 1. Load genome
  
  cat("[2] Loading genome...")
  prepGenome(genome.name, genome.path, genome.prefix, genome.suffix, genome, kmertone.env)
  cat("\n")

  # ---------------- C. GENOMIC COORDINATE --------------------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)

  cat("[3] Loading genomic coordinate table...\n")
  prepGenCoordinate("genomic.coordinate", strand.mode, genome, kmertone.env)
  gc()
  cat("\n")
  
  # add column sequence if strand sensitive
  if (strand.mode == "sensitive") {
    cat("[4] Filtering genomic coordinate table...\n")
    addColumnSequence("genomic.coordinate", genome, kmertone.env)
    filterTable("genomic.coordinate", DNA.pattern, strand.mode, output, kmertone.env) # memory spike here
    fwrite(genomic.coordinate, paste0(output, "/filtered_table.csv"))
  }
  
  assign("genomic.coordinate", genomic.coordinate, envir = globalenv())
  
  # Backup original coordinates
  genomic.coordinate[, c("original_start", "original_end") := list(start, end)]
  
  # ---------------- 4. PRE-ANALYSIS --------------------------------------------------------
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
                      strand.mode, env=kmertone.env)
    
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
  
  # ---------------- 5. KMER EXTRACTION ------------------------------------------------------
  
  cat("[5] Getting case kmers...\n\n")
  kmers <- getCaseKmers("genomic.coordinate", genome, k, DNA.pattern, strand.mode,
               remove.overlaps=TRUE, kmertone.env)
  cat("\n")

  if (!is.null(control.relative.position)) {
    cat("[6] Getting control kmers...\n\n")
    kmers <- getControlKmers("genomic.coordinate", genome, k, DNA.pattern, strand.mode,
                             "kmers", kmertone.env)
  }

  # # ---------------- 6. UPDATE PRE-ANALYSIS ---------------------------------------------------
  # 
  # #GCcontent(dts[1], genome, filename = "GC_after")
  # #Gcontent(dts[1], genome, filename = "G_after")
  # #seqlogo(dts[1], genome, filename = "seqlogo_after")
  # 
  # # ---------------- 7. SCORE -----------------------------------------------------------------

  cat("\n[7] Calculating z score...")
  zScore("kmers", kmertone.env)
  #pValue("kmers") # TBD
  cat("DONE!\n")
 
  # ---------------- THE END ---------------------------------------------------------------
  
  # if (ncpu > 1) {
  #   stopCluster(cl)
  # }
  
  # Save kmers
  if(is.null(DNA.pattern)){
    fwrite(kmers, paste0(output, "/kmers_k-", k, ".csv"))
  } else {
    fwrite(kmers, paste0(output, "/kmers_k-", k, "_pattern-", DNA.pattern, ".csv"))
  }
   
  return(kmers) # kmers table: kmer, count, fold_change, p, z
}