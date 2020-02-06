prepGenome <- function(env) {
  
  # Dependencies:
  #       Kmertone variables: genome.name, [if user's genome: genome.path, genome.suffix]
  

  if (!is.null(genome.path)) {
    
    env$chromosome.names <- list.files(genome.path)[grep(paste0("\\", env$genome.suffix, "$"),
                                                                  list.files(env$genome.path))]
    env$chromosome.names <- gsub(paste0('\\', env$genome.suffix, '$'), '', env$chromosome.names)
    env$genome.name <- "genome"
    
    cat("Genome path is", env$genome.path, "\n")
    cat("Detected chromosome names are\n")
    cat(env$chromosome.names, "\n")
    
  } else if (env$genome.name %in% c("GRCh37", "GRCh38")) {
    
    env$genome.path <- paste0("data/", env$genome.name, "/")
    env$genome.suffix <- ".fa.gz"
    
    env$chromosome.names <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    
  } else {
    
    stop(paste0("Genome ", env$genome.name, " is not available!"))
    
  }
  
  for (chr in env$chromosome.names) {
    
    loadGenome("genome", chr, env$genome.path, env$genome.prefix, env$genome.suffix, env = env, form = "string")
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


