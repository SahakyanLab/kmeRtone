prepGenome <- function(env) {
  
  # Dependencies:
  #       Kmertone variables: genome.name, 
  #                           if user's own genome: genome.path, genome.suffix, genome.prefix OR full.path
  

  if (!is.null(genome.path)) {
    
    env$chromosome.names <- list.files(genome.path)[grep(paste0("\\", env$genome.suffix, "$"),
                                                                  list.files(env$genome.path))]
    env$chromosome.names <- gsub(paste0('\\', env$genome.suffix, '$'), '', env$chromosome.names)
    env$genome.name <- "genome"
    
    cat("\n\nGenome path is", env$genome.path, "\n")
    cat("\nDetected chromosome names are\n")
    cat(env$chromosome.names, "\n")
    
    for (chr in env$chromosome.names) {
      loadGenome("genome", chr, env$genome.path, env$genome.prefix, env$genome.suffix, env = env, form = "string")
    }
    
  } else if (env$genome.name %in% c("GRCh37", "GRCh38", "hg19", "hg38")) {
    
    if (env$genome.name == "hg19") env$genome.name <- "GRCh37"
    if (env$genome.name == "hg38") env$genome.name <- "GRCh38"
    
    genome.path <- paste0("data/genome/", env$genome.name, ".rds")
    env$genome <- readRDS(genome.path)
    env$chromosome.names <- names(env$genome)
    
    cat("DONE!\n")
    
    # add print function
    env$print.genome <- function(obj) print(attr(obj, "length"))
  } else {
    
    stop(paste0("\nGenome ", env$genome.name, " is not available!"))
    
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


