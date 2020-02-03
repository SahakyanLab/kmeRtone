kmertone <- function(genomic.coordinate, genome.name="GRCh37", genome.path=NULL, genome.prefix="", genome.suffix=".fa.gz",
                     DNA.pattern="TT", k=10, expansion.factor=0, control.region,
                     directionality.mode="sensitive", ncpu=1) {
  
  # this is the only function user can call. The rest are internal functions
  # genomic.coordinate     <data.table>     A data.table of genomic coordinate. It can be a list of
  #                                         data.table for replicates. It must contains column 
  #                                         chromosome, start, end, and strand. For "insensitive" mode
  #                                         strand information is not required.
  # genome.name            <string>         Name of available genome: "GRCh37" or "GRCh38"
  # genome.path            <string>         If provided genome is not available, user can input a path
  #                                         to a folder containing chromosome fasta files. The fasta
  #                                         files should be named as chromosome names just like in the
  #                                         genomic coordinate table and has a uniform extension name.
  #                                         The extension can  be specified in flag suffix below. The
  #                                         fasta files are expected to contain only one-line header
  #                                         followed by same-length sequence in every line.
  #                                         e.g. chr1.fasta, chr2.fasta, etc. The 
  # genome.suffix          <string>         An extension name of the fasta files. Compressed file is
  #                                         supported.
  #                                         e.g. .fasta, .fa, .fa.gz, etc.
  # DNA.pattern            <string>         A single DNA pattern. e.g. "G", "TT", "TC", etc. It can be
  #                                         set as NULL if no pattern is desired.
  # k.size                 <integer>        The size of a kmer.
  # control.region         <vector>         A coordinate in a vector format i.e. c(start, end) pointing
  #                                         to where the control kmers should be extracted from. The
  #                                         coordinate is relative to the damage site (upstream and
  #                                         downstream)
  # directionality.mode    <string>         Mode of strand directionality. "sensitive" or "insensitive".
  #                                         In sensitive mode, strand directionality is important. The
  #                                         control kmers are be extracted from the sense strand, i.e.
  #                                         the same strand as the damage. In contrast, in "insensitive"
  #                                         mode, the control kmers are extracted from both sense and
  #                                         antisense strands.
  # ncpu                   <integer>        The number of cpu to use. The default is one.
  
  kmertone.env = environment()
  
  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant functions #########################################################
  source("lib/loadGenomeV2.R")
  source("lib/reverseComplement.R")
  source("lib/distributeChunk.R")
  source("lib/countSlidingBool.R")
  source("lib/scaleGenCoordinate.R")
  
  # Dependant functions from the TrantorR library
  
  
  # Task specific dependant functions
  source("lib/prepGenome.R")
  source("lib/prepGenCoordinate.R")
  source("lib/addColumnSequence.R")
  source("lib/filterTable.R")
  source("lib/damageDistribution.R")
  
  source("lib/checkCoordinate.R")
  source("lib/GCcontent.R")
  source("lib/Gcontent.R")
  source("lib/seqlogo.R")
  source("lib/extractSensitiveKmer.R")
  source("lib/extractInsensitiveKmer.R")
  
  ## Dependant libraries #########################################################
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(stringi)    )

  ## Parallel setup ##############################################################
  if (ncpu > 1) {
    
    suppressPackageStartupMessages( library(foreach)    )
    suppressPackageStartupMessages( library(doParallel) )
    
    cl <- makeCluster(ncpu)
    registerDoParallel(cl)
    
    clusterExport(cl, c("loadGenome", "reverseComplement", "genome.path"))
  }
  
  ## Directory setup #############################################################
  suppressWarnings(dir.create("data"))
  suppressWarnings(dir.create("data/GC"))
  
  
  
  # ---------------- GENOME -------------------------------------------------------------
  
  prepGenome()
  
  # ---------------- GENOMIC COORDINATE --------------------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)

  prepGenCoordinate()
  
  # ---------------- PRE-ANALYSIS --------------------------------------------------------
  # Replicate summary
  # GC and G content at various width
  # Seqlogos
  
  # add column sequence
  addColumnSequence()
  
  # filter table - remove chrM
  damageDistribution(plot=FALSE) # Need more work!
  filterTable()
  
  # damage distribution
  damageDistribution(plot=TRUE)
  plotDistribution() # Need more work!
  
  # G|C and G content
  Ncontent(count.genome = T, count.damage = "retrieve", N="GC") # Need checking
  Ncontent(count.genome = T, count.damage = "retrieve", N="G")  # Need checking
  
  # Plot G|C and G density
  plotDensity() # TBD
  
  # seqlogo
  drawSeqlogo() # TBD
  
  # volcano plot
  plotVolcano() # TBD
  
  # ---------------- KMER EXTRACTION ------------------------------------------------------
  
  if (directionality.mode == "sensitive") {
    
    kmers <<- extractSensitiveKmer() # TBD
    
  } else if (directionality.mode == "insensitive") {
    
    kmers <<- extractInsensitiveKmer() # TBD
    
  } else {
    stop('Please choose either "sensitive" or "insensitive" mode.')
  }
  
  # ---------------- UPDATE PRE-ANALYSIS ---------------------------------------------------
  
  GCcontent(dts[1], genome, filename = "GC_after")
  Gcontent(dts[1], genome, filename = "G_after")
  seqlogo(dts[1], genome, filename = "seqlogo_after")
  
  # ---------------- SCORE -----------------------------------------------------------------
  
  kmers <- p.value(dts[2]) # TBD
  kmers <- z.score(kmers) # TBD
 
  # ---------------- THE END ---------------------------------------------------------------
  
  if (ncpu > 1) {
    stopCluster(cl)
  }
   
  return(kmers) # kmers table: kmer, count, fold_change, p_value, z_score
}