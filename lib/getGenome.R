getGenome <- function(genome.name, genome.path=NULL) {
  
  if (!is.null(genome.path)) {
    
    # convert user input genome fasta file to individual chromosome
    buildGenome()
    
  } else {
    
    # genome in our storage
    genome.list = list(hg19 = hg19.path,
                       hg38 = hg38.path)
    
    genome.path = genome.list[ names(genome.list) %in% genome.name ][[1]]
    
    chromosome.names = c(paste0("chr", 1:23), "chrX", "chrY")
    
    genome = lapply(chromosome.names, function(chr) loadGenome(paste0(genome.path, chr, ".fasta")))
  }
  
  
  return(genome)
}