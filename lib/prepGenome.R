prepGenome <- function() {
  # Global variables: genome.name, [if user's genome: genome.path, genome.suffix]
  
  if (!is.null(genome.path)) {
    
    genome.path <<- genome.path
    chromosome.names <<- list.files(genome.path)[grep(paste0("\\", genome.suffix, "$"), list.files(genome.path))]
    chromosome.names <<- gsub(paste0('\\', genome.suffix, '$'), '', chromosome.names)
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
  
  for (chr in chromosome.names) {
    loadGenome(chr, genome.path, genome.prefix, genome.suffix, kmertone.env, form = "string")
  }
}