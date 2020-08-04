calculateOptimumK <- function(genomic.coordinate, genome, set.k, DNA.pattern,
                              strand.mode, env=parent.frame()) {
  # To calculate optimum k value. The calculation is based on
  # Claudia's thesis equation 3 and 4.
  # This is a stand-alone function - separate from kmertone function.
  
  # Dependencies:
  #     Packages: data.table, stringi
  #     Functions: prepGenome, prepGenCoordinate, addColumnSequence, filterTable, getCaseKmers
  #     Information: 1) proportion of precalculated A, C, G, and T of the whole genome (hg19)
  
  set.k <- sort(set.k)
  base.k <- if(k[1]%%2 == 0) k[1]-2 else if(k[1]%%2 > 0) k[1]-2
  if(base.k > 0) set.k <- c(base.k, set.k)
  
  ## Dependant libraries #########################################################
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(  stringi ) )
  
  ## Dependant functions #########################################################
  source("lib/loadGenomeV2.R")
  source("lib/reverseComplement.R")
  source("lib/distributeChunk2.R")
  source("lib/scaleGenCoordinate.R")
  source("lib/trimGenCoordinates.R")
  source("lib/expandGenCoordinate.R")
  source("lib/extractKmers.R")
  source("lib/countReverseComplementKmers.R")
  source("lib/removeAllOverlaps.R")
  
  # Task specific dependant functions
  source("lib/02_prepGenome.R", local = TRUE)
  source("lib/03a_prepGenCoordinate.R", local = TRUE)
  source("lib/03b_addColumnSequence.R", local = TRUE)
  source("lib/03c_filterTable.R", local = TRUE)
  source("lib/07_getCaseKmers.R", local = TRUE)
  
  main.env <- environment()

  ####### PREPARATION ###################################################################
  # ---------------- GENOME -------------------------------------------------------------
  # 1. Load genome
  
  cat("[1] Loading genome...")
  prepGenome(genome.name, genome.path, genome.prefix, genome.suffix, genome, main.env)
  cat("\n")
  
  # ---------------- GENOMIC COORDINATE --------------------------------------------------
  # 1. Rename columns
  # 2. Combine replicates (if any)
  
  cat("[2] Loading genomic coordinate table...\n")
  prepGenCoordinate("genomic.coordinate", strand.mode, genome, main.env)
  #fwrite(genomic.coordinate, "data/merged_table.csv")
  gc()
  cat("\n")
  
  # add column sequence if strand sensitive
  if (strand.mode == "sensitive") {
    cat("[3] Filtering genomic coordinate table...\n")
    addColumnSequence("genomic.coordinate", genome, main.env)
    filterTable("genomic.coordinate", DNA.pattern, strand.mode, main.env) # memory spike here
    #fwrite(genomic.coordinate, "data/filtered_table.csv")
  }
  
  # Backup original coordinates
  genomic.coordinate[, c("original_start", "original_end") := list(start, end)]
  
  # ---------------- MAIN FUNCTION --------------------------------------------------------
  
  cat("[4] Extracting case kmers...\n")
  kmers.list <- lapply(set.k, function(k) {
    
    cat("k =", k, "\n")
    kmers <- getCaseKmers("genomic.coordinate", genome, k, DNA.pattern, strand.mode, 
                          remove.overlaps=FALSE, env=main.env)
    main.env[["genomic.coordinate"]][, c("start", "end") := list(original_start, original_end)]
    
    return(kmers)
  })
  
  names(kmers.list) <- set.k
  
  # for (k in set.k) {
  #   fwrite(kmers.list[[as.character(k)]], paste0("data/case-kmers_k-", k, ".csv"))
  # }
  
  cat("[5] Calculating optimal k...\n")
  
  q <- sapply(set.k[-1], function(k) {
    
    q <- kmers.list[[as.character(k)]][, {
      
      total.case <- sum(count)
      
      nt.init <- stri_sub(kmer, 1, 1)
      seq.mid <- stri_sub(kmer, 2, k-1)
      nt.last <- stri_sub(kmer, k, k)
      
      # take proportion from previous k
      p.mid <- kmers.list[[as.character(k-2)]][, {
        
        names(count) <- kmer
        count[seq.mid]/sum(count)
        
        }]
      
      # proportion of precalculated A, C, G, and T of the whole genome (hg19)
      p.nt <- c(A=0.295, T=0.295, C=0.205, G=0.205)
      
      p.init <- p.nt[nt.init]
      p.last <- p.nt[nt.last]
      
      p = p.init * p.mid * p.last
      
      # equation 3 - to compute chi square
      Q = total.case * sum ( (count/total.case - p)**2 / p )
      
      q = ( Q - (4**(k-2) - 1) ) / sqrt( 2 * (4**(k-2) - 1) )
      
      q
    }]
    
    names(q) <- as.character(k)
    
    return(q)
  })
  
}