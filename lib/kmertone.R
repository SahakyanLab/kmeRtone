kmertone <- function(genomic.coordinate, genome.name, genome.path=NULL, damage.pattern,
                     k.size, control.region, directionality.mode, ncpu=1) {
  
  # this is the only function user can call. The rest are internal functions
  # genomic.coordinate     A data.table of genomic coordinate
  # genome.name            Name of available genome: hg19 or hg38
  # genome.path            If provided genome is not available, user can input a single
  #                        fasta file containing all chromosome sequences with chromosome
  #                        names as the headers.
  # damage.pattern         A single damage pattern in a string format.
  # k.size                 The size of a kmer
  # control.region         A coordinate in a vector format i.e. c(start, end) pointing to
  #                        where control kmers should be extracted from.
  # directionality.mode    Should the analysis be in a "sensitive" or "insensitive" mode for
  #                        strand directionalty?
  
  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant functions #########################################################
  source("lib/fsa2fastas.R")
  source("lib/loadGenomeAsList.R")
  source("lib/genomeSequenceMapping.R")
  source("lib/reverseComplement.R")
  
  # Dependant functions from the TrantorR library
  source(paste0(TrantoRLib, "GEN_readfasta.R"))
  source(paste0(TrantoRLib, "UTIL_readLinesFast.R"))
  source(paste0(TrantoRLib, "GEN_loadGenome.R"))
  
  # Task specific dependant functions
  source("lib/getGenome.R")
  source("lib/checkCoordinate.R")
  source("lib/GCcontent.R")
  source("lib/Gcontent.R")
  source("lib/seqlogo.R")
  source("lib/extractSensitiveKmer.R")
  source("lib/extractInsensitiveKmer.R")
  
  ## Dependant libraries #########################################################
  suppressPackageStartupMessages( library(data.table) )
  
  ## Parallel setup ##############################################################
  if (ncpu > 1) {
    source("lib/guf_doParallelsetup_nofun.R")
  }
  
  
  ## Directory setup #############################################################
  suppressWarnings(dir.create("data"))
  
  # ---------------- GENOME -------------------------------------------------------------
  
  genome = getGenome(genome.name, genome.path)
  
  # ---------------- GENOMIC COORDINATE --------------------------------------------------

  # use dt for variable name and enforce to use data.table
  dt = as.data.table(genomic.coordinate)
  
  # check the table
  checkCoordinate(dt, genome, damage.pattern)
  
  # ---------------- PRE-ANALYSIS --------------------------------------------------------
  # GC and G content at various width
  # Seqlogos
  
  GCcontent(dt, genome, filename = "GC_before")
  Gcontent(dt, genome, filename = "G_before")
  seqlogo(dt, genome, filename = "seqlogo_before")
  
  # ---------------- KMER EXTRACTION ------------------------------------------------------
  
  if (directionality.mode == "sensitive") {
    
    dts <- extractSensitiveKmer(dt, genome, damage.pattern, k.size, control.region)
    
  } else if (directionality.mode == "insensitive") {
    
    dts <- extractInsensitiveKmer(dt, genome, damage.pattern, k.size, control.region)
    
  } else {
    stop('Please only choose "sensitive" or "insensitive" mode')
  }
  
  # ---------------- UPDATE PRE-ANALYSIS ---------------------------------------------------
  
  GCcontent(dts[1], genome, filename = "GC_after")
  Gcontent(dts[1], genome, filename = "G_after")
  seqlogo(dts[1], genome, filename = "seqlogo_after")
  
  # ---------------- SCORE -----------------------------------------------------------------
  
  kmers <- p.value(dts[2])
  kmers <- z.score(kmers)
  
  # ---------------- PLOT ------------------------------------------------------------------
  
  volcanoPlot(kmers)
  
  return(kmers) # kmers table: kmer, count, fold_change, p_value, z_score
}