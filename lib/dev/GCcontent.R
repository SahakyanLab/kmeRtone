GCcontent <- function(chromosome.names, genome.name, count.genome=F, count.damage=F) {
  # To calculate GC percentage of damage and damage pattern in the genome
  #     at several window size.
  # 1) To load pre-calculated GC percentage of the genome or calculate from scratch
  # 3) To plot the density of GC distribution.
  
  chromosome.names <- c(paste0("chr", 1:22), "chrX", "chrY")
  
  width <- c(100, 500, 1000, 5000, 10000, 50000, 100000)
  chromosomes.size <- sapply(chromosome.names, function(chr) readGenome(chr, size = TRUE))
  
  # -------------------------------------------------------------------------------------------
  # GC content of genome
  if (count.genome == TRUE) {
    
    for (chr in chromosome.names) {
      
      # get chromosome sequence
      chromosome.seq <- readGenome(chr, 1, chromosomes.size[[chr]])
      
      # get GC boolean
      GC.bool <- chromosome.seq == "G" | chromosome.seq == "C"
      
      # we don't need the sequence anymore. discard!
      rm(chromosome.seq)
      
      cat("\nCounting GC of", chr, "at", as.integer(width[1]), "width\n")
      GC.count <- countSlidingBool(GC.bool, width[1], NCPU)
      
      # save result
      fwrite(list(GC.count), paste0("data/GC/", genome.name, "_", chr, "_GC_width_", as.integer(width[1]), ".csv"))
      
      GC.count.1 <- GC.count
      w1 <- width[1]
      
      for (w in width[-1]) {
        cat("Counting GC of", chr, "at", as.integer(w), "width\n")
        
        GC.count <- scaleCountSliding(GC.count.1, w1, w)
        
        GC.count.1 <- GC.count
        w1 <- w
        
        # save result
        fwrite(list(GC.count), paste0("data/GC/", genome.name, "_", chr, "_GC_width_", as.integer(w), ".csv"))
      }
    }
    rm(GC.count, GC.bool, GC.count.1)
  }
  
  # -------------------------------------------------------------------------------------------
  # GC content of the damage coordinate
  # Retrieve the pre-calculated GC percentage of the genome or recalculate from scratch
  
  if (count.damage == FALSE) {
    
    for (w in width) {
      
      dt.w <- scaleGenomicCoordinate(dt[, 1:4], scale = w, side = both)
      setkey(dt.w, chromosome, start, end)
      
      dt.w[, {
        
        # call a function to retrieve GC
        GC.count <- retrieveGC(chromosome, start-w, end+w, w)
        
        # save result
        fwrite(list(GC.count), paste0("data/GC/damage_", chromosome, "_GC_width_", w, ".csv"))
        
      }, by = chromosome] 
    }
    
  } else if (count.damage == TRUE) {
    
    for (w in width[-1]) {
      
      dt.w <- scaleGenomicCoordinate(dt[, 1:4], scale = w, side = both)
      setkey(dt.w, chromosome, start, end)
      
      dt.w[, {
        
        # get dna sequence
        dna.seqs <- readGenome(chromosome, start, end)
        dna.seq <- unlist(dna.seqs)
        dna.seq <- matrix(dna.seq, ncol=.N)
        
        # count GC
        GC.count <- colSums(dna.seqs == "G" | dna.seqs == "C")
        
        fwrite(list(GC.count), paste0("data/GC/damage_", chromosome, "_GC_width_", w, ".csv"))
        
      }, by = chromosome]
      
    }
  }
  
  # -------------------------------------------------------------------------------------------
  # Plot GC density
  
  chromosome.names <- dt[, unique(chromosome)]
  
  for (w in width) {
    
    for (chr in chromosome.names) {
      
      # load genome GC content
      genome.GC <- fread("data/GC/", genome.name, "_GC_width_", w, ".csv.gz")[[1]]
      damage.GC <- fread("data/GC/damage", "_GC_width_", w, ".csv.gz")[[1]]
      
      gen.density <- density(genome.GC)
      dmg.density <- density(damage.GC)
      
      plot(gen.density)
      
      line(dmg.density)
      
     # save plot 
      
    }
  }
}






