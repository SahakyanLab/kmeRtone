kmertone <- function(genomic.coordinate, genome.name="GRCh37", strand.mode,
                     k, control.relative.position, DNA.pattern=NULL, ncpu=1,
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
  # DNA.pattern                  <string>         A single DNA pattern. e.g. "G", "TT", "TC", etc. It can be
  #                                               set as NULL if no pattern is desired.
  # k.size                       <numeric>        The size of a kmer.
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
  source("lib/countSlidingBool.R")
  source("lib/scaleGenCoordinate.R")
  source("lib/trimGenCoordinates.R")
  source("lib/mergeGenCoordinate.R")
  source("lib/expandGenCoordinate.R")
  source("lib/extractKmers.R")
  
  # Dependant functions from the TrantorR library
  source("lib/TrantoRext/GEN_getSeqLogo.R", local = TRUE)
  
  # Task specific dependant functions
  source("lib/inputChecking.R", local = TRUE)
  source("lib/prepGenome.R", local = TRUE)
  source("lib/prepGenCoordinate.R", local = TRUE)
  source("lib/addColumnSequence.R", local = TRUE)
  source("lib/filterTable.R", local = TRUE)
  source("lib/damageDistribution.R", local = TRUE)
  source("lib/Ncontent.R", local = TRUE)
  source("lib/removeCaseZone.R")
  source("lib/getKmers.R", local = TRUE)
  source("lib/zScore.R", local = TRUE)
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
  suppressWarnings(dir.create("data/GC"))
  
  
  # ---------------- INPUT CHECKING -----------------------------------------------------
  
  cat("[1] Checking inputs...")
  inputChecking(kmertone.env)
  cat("DONE!\n")
  
  # ---------------- GENOME -------------------------------------------------------------
  # 1. Load genome
  
  cat("[2] Loading genome...")
  prepGenome(kmertone.env)
  cat("DONE!\n")
  
  # ---------------- GENOMIC COORDINATE --------------------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)

  cat("[3] Loading genomic coordinate table...")
  prepGenCoordinate(kmertone.env)
  cat("DONE!\n")
  
  # ---------------- PRE-ANALYSIS --------------------------------------------------------
  # Replicate summary
  # GC and G content at various width
  # Seqlogos
  
  # add column sequence if strand sensitive
  if (strand.mode == "sensitive") {
    cat("[4] Filtering genomic coordinate table...\n")
    addColumnSequence(kmertone.env)
    filterTable(kmertone.env)
  }
  
  # filter table - remove chrM
  #damageDistribution() # Need more work!
  
  
  # damage distribution
  #damageDistribution(plot=TRUE)
  #plotDistribution() # Need more work!
  
  # G|C and G content
  #Ncontent(count.genome = T, count.damage = "retrieve", N="GC") # Need checking
  #Ncontent(count.genome = T, count.damage = "retrieve", N="G")  # Need checking
  
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
  
  cat("[6] Calculating z score...")
  zScore("kmers")
  #pValue("kmers") # TBD
  cat("DONE!\n")
 
  # ---------------- THE END ---------------------------------------------------------------
  
  if (ncpu > 1) {
    stopCluster(cl)
  }
   
  return(kmers) # kmers table: kmer, count, fold_change, p_value, z_score
}