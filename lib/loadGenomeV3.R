loadGenome <- function(genome.path, genome.name, env=parent.frame(), save.as.RDS=FALSE){
  # Load and build a genome-class object. I had tested S3 class vs S4 class (DNAStringSet) in term of speed.
  # S3 class is significantly faster for subsetting.
  
  # genome.path     <string>    A path to a genome folder containing chromosome single-header fasta files.
  #                             Each header should contain appropriate chromosome name.
  #                             The fasta file should follow UCSC Genome convention.
  
  # Dependency
  #    Packages: stringi, data.table,
  
  # List chromosome fasta files
  chr.filenames <- list.files(genome.path, "(chr[0-9A-Za-z_]+|mito.+)\\.(fa|fna|fasta)")

  # Load all genome
  genome <- sapply(chr.filenames, function(chr.filename){
    
    chr.path <- paste0(genome.path, "/", chr.filename)
    chr <- readLines(chr.path)
    chr.name <- sub("^>", "", chr[1])
    chr.seq <- stri_c(chr[-1], collapse = "")
    chr.seq <- stri_trans_toupper(chr.seq)
    names(chr.seq) <- chr.name

    return(chr.seq)
  }, USE.NAMES = FALSE)
  
  genome <- as.list(genome)
  
  # Assign class genome
  class(genome) <- "genome"
  env$print.genome <- function(obj) print(attr(obj, "length"))
  
  # Add attribute length
  attr(genome, "length") <- sapply(genome, nchar)
  
  # Save as rds file
  if(save.as.RDS){
    saveRDS(genome, paste0(genome.path, "/", genome.name, ".rds"))
  }
  
  return(genome)
}