prepGenome <- function(genome.name, genome.path, genome.prefix,
                       genome.suffix, genome, env) {
  
  # Dependencies:
  #       Kmertone variables: genome.name, 
  #                           if user's own genome: genome.path, genome.suffix, genome.prefix OR full.path
  
  if (class(genome) == "genome") {
    
    # add print function
    env$print.genome <- function(obj) print(attr(obj, "length"))
    env$genome <- genome
    
    cat("Genome is already loaded.\n\n")
    print(genome)

  } else if (!is.null(genome.path)) {
    
    chromosome.names <- list.files(genome.path)[grep(paste0("\\", genome.suffix, "$"),
                                                                  list.files(genome.path))]
    chromosome.names <- gsub(paste0('\\', genome.suffix, '$'), '', chromosome.names)
    genome.name <- "genome"
    
    cat("\n\nGenome path is", genome.path, "\n")
    cat("\nDetected chromosome names are\n")
    cat(chromosome.names, "\n")
    
    for (chr in chromosome.names) {
      loadGenome("genome", chr, genome.path, genome.prefix, genome.suffix, env = env, form = "string")
    }
    
  } else if (genome.name %in% c("GRCh37", "GRCh38", "hg19", "hg38")) {
    
    if (genome.name == "hg19") genome.name <- "GRCh37"
    if (genome.name == "hg38") genome.name <- "GRCh38"
    
    genome.path <- paste0("data/genome/", genome.name, ".rds")
    env$genome <- readRDS(genome.path)
    chromosome.names <- names(env$genome)
    
    cat("DONE!\n")
    
    # add print function
    env$print.genome <- function(obj) print(attr(obj, "length"))
  } else {
    
    stop(paste0("\nGenome ", genome.name, " is not available!"))
    
  }
}


# # build a S3-class genome containing reference class chromosome
# 
# setChromosome <- setRefClass("chromosome", fields = list(sequence = "character"))
# 
# setChromosome$methods(len = function() length <<- nchar(sequence))
# 
# # printing function
# setChromosome$methods(show = function(obj) {
#   cat("Sequence:",
#       substr(sequence, 1, 3), "...",
#       substr(sequence, nchar(sequence)-10, nchar(sequence)), "\n")
#   cat("Length:", length, "\n")
# })
# 
# chr1 <- setChromosome(sequence = "AAAAAAAAAAAAA", length = nchar(sequence))


