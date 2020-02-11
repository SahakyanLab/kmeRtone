kmertone <- function(genomic.coordinate, genome.name="GRCh37", strand.mode,
                     k, control.relative.position, DNA.pattern=NULL, ncpu=1,
                     genomic.coordinate.background = NULL,
                     genome.path=NULL, genome.prefix="", genome.suffix=NULL) {
  
  # this is the only function user can call. The rest are internal functions
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
  # control.relative.position    <vector>         A coordinate in a vector format i.e. c(start, end) pointing
  #                                               to where the control kmers should be extracted from. The
  #                                               coordinate is relative to the damage site (upstream and
  #                                               downstream)
  # strand.mode                  <string>         Mode of strand: "sensitive" or "insensitive".
  #                                               In "sensitive" mode, strand information is important. The
  #                                               kmers are extracted from the sense strand only.
  #                                               In contrast, in "insensitive" mode, the kmers are extracted
  #                                               from both plus and minus strands.
  # ncpu                         <numeric>        The number of cpu to use. The default is one.
  
  # library(data.table)
  # genomic.coordinate <- fread("data/table.csv")
  # genome.name="GRCh37";genome.path=NULL;genome.prefix=""; genome.suffix=".fa.gz"
  # DNA.pattern="TT"; k=10; control.relative.position=c(80,500)
  # strand.mode="sensitive"; ncpu=1

  kmertone.env = environment()
  
  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant functions #########################################################
  source("lib/loadGenomeV2.R")
  source("lib/reverseComplement.R")
  source("lib/distributeChunk.R")
  source("lib/distributeChunk2.R")
  source("lib/countSlidingBool.R")
  source("lib/scaleGenCoordinate.R")
  source("lib/trimGenCoordinates.R")
  source("lib/mergeGenCoordinate.R")
  source("lib/expandGenCoordinate.R")
  source("lib/extractKmers.R")
  source("lib/removeAllOverlaps.R")
  source("lib/removeCaseZone.R")
  
  # Dependant functions from the TrantorR library
  source("lib/TrantoRext/GEN_getSeqLogo.R", local = TRUE)
  
  # Task specific dependant functions
  source("lib/01_inputChecking.R", local = TRUE)
  source("lib/02_prepGenome.R", local = TRUE)
  source("lib/03a_prepGenCoordinate.R", local = TRUE)
  source("lib/03b_addColumnSequence.R", local = TRUE)
  source("lib/03c_filterTable.R", local = TRUE)
  source("lib/04_damageDistribution.R", local = TRUE)
  source("lib/05_GCcontent.R", local = TRUE)
  source("lib/07_getKmers.R", local = TRUE)
  source("lib/08_zScore.R", local = TRUE)
  #source("lib/seqlogo.R", local = TRUE)

  ## Dependant libraries #########################################################
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(stringi)    )

  ## Parallel setup ##############################################################
  if (ncpu > 1) {
    
    suppressPackageStartupMessages( library(foreach)    )
    suppressPackageStartupMessages( library(doParallel) )
    
    cl <- makeCluster(ncpu)
    registerDoParallel(cl)

  }
  
  ## Directory setup #############################################################
  suppressWarnings(dir.create("data"))

  
  # ---------------- INPUT CHECKING -----------------------------------------------------
  
  cat("[1] Checking inputs...")
  inputChecking(kmertone.env)
  cat("DONE!\n")
  
  # ---------------- GENOME -------------------------------------------------------------
  # 1. Load genome
  
  cat("[2] Loading genome...")
  prepGenome(kmertone.env)
  
  # ---------------- GENOMIC COORDINATE --------------------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)

  cat("[3] Loading genomic coordinate table...\n")
  prepGenCoordinate(kmertone.env)
  #fwrite(genomic.coordinate, "data/merged_table.csv")
  gc()
  
  # add column sequence if strand sensitive
  if (strand.mode == "sensitive") {
    cat("[4] Filtering genomic coordinate table...\n")
    addColumnSequence(kmertone.env)
    gc()
    filterTable(kmertone.env) # memory spike here
    #fwrite(genomic.coordinate, "data/filtered_table.csv")
  }
  
  # ---------------- PRE-ANALYSIS --------------------------------------------------------
  # Replicate summary
  # GC and G content at various width
  # Seqlogos
  
  # damage distribution
  #damageDistribution(plot=TRUE)
  #plotDistribution() # Need more work!
  
  # G|C and G content
  GCcontent(env)

  # Plot G|C and G density
  #plotDensity() # TBD
  
  # seqlogo
  #drawSeqlogo() # TBD
  
  # volcano plot
  #plotVolcano() # TBD
  
  # ---------------- KMER EXTRACTION ------------------------------------------------------
  
  cat("[5] Getting kmers...\n\n")
  getKmers(kmertone.env)
  
  # ---------------- UPDATE PRE-ANALYSIS ---------------------------------------------------
  
  #GCcontent(dts[1], genome, filename = "GC_after")
  #Gcontent(dts[1], genome, filename = "G_after")
  #seqlogo(dts[1], genome, filename = "seqlogo_after")
  
  # ---------------- SCORE -----------------------------------------------------------------
  
  cat("\n[6] Calculating z score...")
  zScore("kmers")
  #pValue("kmers") # TBD
  cat("DONE!\n")
 
  # ---------------- THE END ---------------------------------------------------------------
  
  if (ncpu > 1) {
    stopCluster(cl)
  }
   
  return(kmers) # kmers table: kmer, count, fold_change, p, z
}