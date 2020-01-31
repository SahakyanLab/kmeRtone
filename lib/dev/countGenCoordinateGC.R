countGenCoordinateGC <- function(genomic.coordinate, genome.path, ncpu=1) {
  # This function calculate percentage of GC of a given genomic coordinate table
  
  # Dependencies: data.table; readGenome
  
  dt <- genomic.coordinate
  setkey(dt, chromosome, start, end)
  
  if (ncpu == 1) {
    
    dt[, GC := {
      
      dna.seq <- readGenome(chromosome, min(start), max(end))
      
      GC.bool <- dna.seq == "G" | dna.seq == "C"
      
      rm(dna.seq)
      
      # offset the coordinate
      start.min <- min(start)
      start <- start-start.min+1
      end <- end-start.min+1
      
      GC.percent <- mapply(function(start, end) sum(GC.bool[start:end]) / (end-start+1) * 100,
             start, end)
      
      GC.percent
    }, by = .(chromosome)]
    
  } else {
    
    chunks <- distributeChunk(nrow(dt), ncpu)
    
    dt.chunk <- lapply(1:ncpu, function(i) dt[ chunks$start[i]:chunks$end[i] ])

    GC.percent <- foreach(dt=dt.chunk, .combine = "c", .packages = "data.table", .export = c("readGenome", "genome.PATH"),
                          .noexport = c("dt", "dt.chunk")) %dopar% {
      
        GC.percent <- dt[, GC <- {
        
        dna.seq <- readGenome(chromosome, min(start), max(end))
        
        GC.bool <- dna.seq == "G" | dna.seq == "C"
        
        rm(dna.seq)
        
        # offset the coordinate
        start.min <- min(start)
        start <- start-start.min+1
        end <- end-start.min+1
        
        GC.percent <- mapply(function(start, end) sum(GC.bool[start:end]) / (end-start+1) * 100,
                             start, end)
        
        GC.percent
      }, by = chromosome][[2]]
      
      return(GC.percent)
    }
    
    dt[, GC := GC.percent]
  }
  
  return(dt)
}