kmertone <- function(genomic.coordinate, genome.name="GRCh37", genome.path=NULL, damage.pattern,
                     k.size, control.region, directionality.mode, ncpu=1) {
  
  # this is the only function user can call. The rest are internal functions
  # genomic.coordinate     A data.table of genomic coordinate
  #                        It can be a list of data.table for replicates
  # genome.name            Name of available genome: GRCh37 or GRCh38
  # genome.path            If provided genome is not available, user can input a path to
  #                        a folder containing chromosome fasta files.
  # damage.pattern         A single damage pattern in a string format.
  # k.size                 The size of a kmer
  # control.region         A coordinate in a vector format i.e. c(start, end) pointing to
  #                        where control kmers should be extracted from.
  # directionality.mode    Should the analysis be in a "sensitive" or "insensitive" mode for
  #                        strand directionalty?
  
  # genome.path
  # To make thing simpler, a genome folder must contain fasta files with chromosome name
  # as its filename. The chromosome name must be the same like in the genomic coordinate
  # table. A fasta file must contain a single header with consistent number of bases each 
  # line. It can be in a compressed or non-compressed form.
  
  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant functions #########################################################
  source("lib/readGenome.R")
  source("lib/reverseComplement.R")
  source("lib/addColumnSequence.R")
  Rcpp::sourceCpp("lib/cpp/countSlidingBool.cpp")
  Rcpp::sourceCpp("lib/cpp/scaleCountSliding.cpp")
  source("lib/distributeChunk.R")
  source("lib/countSlidingBool.R")
  source("lib/scaleCountSliding.R")
  
  # Dependant functions from the TrantorR library
  
  
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
  suppressPackageStartupMessages( library(stringi)    )

  ## Parallel setup ##############################################################
  NCPU = ncpu
  if (ncpu > 1) {
    
    suppressPackageStartupMessages( library(foreach)    )
    suppressPackageStartupMessages( library(doParallel) )
    
    cl <- makeCluster(NCPU)
    registerDoParallel(cl)
  }
  
  
  ## Directory setup #############################################################
  suppressWarnings(dir.create("data"))
  suppressWarnings(dir.create("data/GC"))
  
  
  
  # ---------------- GENOME -------------------------------------------------------------
  
  if (genome.name %in% c("GRCh37", "GRCh38")) {
    genome.PATH = paste0("data/", genome.name, "/")
  } else if (!is.null(genome.path)) {
    genome.PATH = genome.path
  } else {
    stop(paste0("Genome", genome.name, "is not available!"))
  }
  
  # only copy necessary objects to cluster
  if (ncpu > 1) {
    clusterExport(cl, c("data.table", "readGenome", "reverseComplement", "genome.PATH"))
  }
  
  # ---------------- GENOMIC COORDINATE --------------------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)

  dt <- genomic.coordinate
  
  # add column sequence
  if (class(dt) == "list") {
    
    cat("Detecting", length(dt), "replicates\n")
    
    # rename columns
    for (dt.rep in dt) {
      colnames(dt.rep) <- c("chromosome", "start", "end", "strand")
    }
    
    # combine replicate
    dt <- combineReplicate(dt)

  } else {
    
    colnames(dt) <- c("chromosome", "start", "end", "strand")
    
  }
  
  chromosome.names <- dt[, unique(chromosome)]
  
  # ---------------- PRE-ANALYSIS --------------------------------------------------------
  # Replicate summary
  # GC and G content at various width
  # Seqlogos
  
  # add column sequence
  addColumnSequence(dt, NCPU)
  
  # filter table
  dt <- filterTable(dt, damage.pattern)
  
  # damage distribution
  damageDistribution(dt)
  
  #
  GCcontent(dt, filename = "GC_before")
  
  
  Gcontent(dt, filename = "G_before")
  seqlogo(dt, filename = "seqlogo_before")
  
  # ---------------- KMER EXTRACTION ------------------------------------------------------
  
  if (directionality.mode == "sensitive") {
    
    kmers <- extractSensitiveKmer(dt, genome, damage.pattern, k.size, control.region)
    
  } else if (directionality.mode == "insensitive") {
    
    kmers <- extractInsensitiveKmer(dt, genome, damage.pattern, k.size, control.region)
    
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
 
  
  # ---------------- THE END ---------------------------------------------------------------
  
  if (ncpu > 1) {
    stopCluster(cl)
  }
   
  return(kmers) # kmers table: kmer, count, fold_change, p_value, z_score
}