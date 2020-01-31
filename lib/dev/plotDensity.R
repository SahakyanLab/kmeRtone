plotDensity <- function() {
  
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