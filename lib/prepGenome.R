prepGenome <- function() {
  # Global variables: genome.name, [if user's genome: genome.path, genome.suffix]
  
  if (!is.null(genome.path)) {
    
    genome.path <<- genome.path
    chromosome.names <<- sub(paste0('\\', genome.sufix, '$'), '', list.files(genome.path))
    genome.name <<- "genome"
    
    cat("Genome path is", genome.path, "\n")
    cat("Detected chromosome names are\n")
    cat(chromosome.names, "\n")
    
  } else if (genome.name %in% c("GRCh37", "GRCh38")) {
    
    genome.path <<- paste0("data/", genome.name, "/")
    
    chromosome.names <<- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    
  } else {
    
    stop(paste0("Genome ", genome.name, " is not available!"))
    
  }
  
}