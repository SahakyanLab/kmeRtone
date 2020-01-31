Ncontent <- function(count.genome=F, count.damage="retrieve", N="N") {
  # To calculate N percentage of damage and damage pattern in the genome
  #     at several window size. N is either G or GC
  # 1) To load pre-calculated N percentage of the genome or calculate from scratch
  # 3) To plot the density of N distribution.
  
  # Dependencies:
  #     Global variables: genomic.coordinate, chromosome.names, genome.path, ncpu
  #     Packages        : data.table, foreach, doParallel
  #     Functions       : readGenome, reverseComplement
  
  if (nchar(N) == 1) {
    N.bool.eval <- paste0("dna.seq == ", N)
  } else {
    base <- strsplit(N, "")[[1]]
    N.bool.eval <- paste0("dna.seq == ", base[1], " | ", "dna.seq == ", base[2])
  }
  

  width <- c(100, 500, 1000, 5000, 10000, 50000, 100000)
  chromosome.sizes <- sapply(chromosome.names, function(chr) readGenome(chr, size = TRUE))
  
  # -------------------------------------------------------------------------------------------
  # N content of genome
  if (count.genome == TRUE) {
    
    for (chr in chromosome.names) {
      
      # get chromosome sequence
      chromosome.seq <- readGenome(chr, 1, chromosomes.size[[chr]])
      
      # get N boolean
      N.bool <- eval(parse(text = N.bool.eval))
      
      # we don't need the sequence anymore. discard!
      rm(chromosome.seq)
      
      cat("\nCounting", N, " of", chr, "at", as.integer(width[1]), "width\n")
      N.count <- countSlidingBool(N.bool, width[1], ncpu)
      
      # save result
      fwrite(list(N.count), sprintf("data/%s/%s_%s_%s_width_%d.csv", N, genome.name, chr, N, as.integer(width[1])))

      N.count.1 <- N.count
      w1 <- width[1]
      
      for (w in width[-1]) {
        cat("Counting N of", chr, "at", as.integer(w), "width\n")
        
        N.count <- scaleCountSliding(N.count.1, w1, w)
        
        N.count.1 <- N.count
        w1 <- w
        
        # save result
        fwrite(list(N.count), sprintf("data/%s/%s_%s_%s_width_%d.csv", N, genome.name, chr, N, as.integer(width[1])))
      }
    }
    rm(N.count, N.bool, N.count.1)
  }
  
  # -------------------------------------------------------------------------------------------
  # N content of the damage coordinate
  # Retrieve the pre-calculated N percentage of the genome or recalculate from scratch
  
  dt <- genomic.coordinate
  setkey(dt, chromosome)
  
  if (count.damage == "retrieve") {
    
    for (w in width) {
      
      dt.w <- scaleGenCoordinate(dt[, 1:4], chromosome.sizes, scale = w/2, side = "both")
      setkey(dt.w, chromosome, start, end)
      
      dt.w[, {
        
        dmg.width <- end[1]-start[1]+1
        
        # call a function to retrieve N count from the genome calculation
        N.count <- retrieveCount(genome.name, chromosome, start, round(dmg.width, -1), "N")
        
        # add extra tail sequence if the width is not exactly the same
        if (w != dmg.width) {
          dna.seq <- matrix(unlist(readGenome(chromosome, start+w, end)), nrow = dmg.width-w)
          N.count <- N.count + colSums(dna.seq == "G" | dna.seq == "C")
        }
        
        # save result
        #fwrite(list(N.count), sprintf("data/%s/%s_%s_%s_width_%d.csv", N, genome.name, chr, N, as.integer(width[1])))
        
      }, by = chromosome] 
    }
    
  } else if (count.damage == "scratch") {
    
    for (w in width) {
      
      # scale up
      dt.w <- scaleGenCoordinate(dt[, 1:4], chromosome.sizes, scale = w/2, side = "both")
      setkey(dt.w, chromosome, start, end)
      
      # can spike up memory depending on the sequence length
      #countGenCoordinateN(dt.w, "")
      
      # if memory exhaustive, use this
      dt.w[, N := {
        
        start.min <- min(start)
        end.max <- max(end)
        dmg.width <- end[1]-start[1]+1
        
        dna.seq <- readGenome(chromosome, start.min, end.max)

        N.bool <- eval(parse(text = N.bool.eval))
        rm(dna.seq)
        
        # calculate sliding window of 100 width
        N.count.w100 <- countSlidingBool(N.bool, 100, ncpu)
        rm(N.bool)

        scale <- w / 100 - 1
        
        # offset the coordinate
        start <- start - start.min + 1

        # assign the 100 parts
        N.count <- mapply(function(start) sum(N.count.w100[seq(start, start+scale*100, 100)]), start)
        
        # reset the offset
        start <- start + start.min - 1
        
        # add extra tail sequence if the width is not exactly the same
        if (w != dmg.width) {
          dna.seq <- matrix(unlist(readGenome(chromosome, start+w, end)), nrow = dmg.width-w)
          N.count <- N.count + colSums(dna.seq == "G" | dna.seq == "C")
        }

        #fwrite(list(N.count), sprintf("data/%s/%s_%s_%s_width_%d.csv", N, genome.name, chr, N, as.integer(width[1])))
        N.count
      }, by = chromosome]
    }
  }
}

