loadGenome <- function(chromosome, genome.path, prefix, suffix,
                       FULL.PATH=NULL, form="vector", case="upper", vector.factor=T) {
  
  
  # Dependencies:
  #      Global variables: genome, if not exist create one in parent frame variable
  
  if (!"genome" %in% ls(envir = parent.frame())) {
    assign("genome", list(), parent.frame())
  }
  
  if (is.null(FULL.PATH)) {
    # get chromosome path
    chromosome.path <- paste0(genome.path, prefix, chromosome, suffix) 
  } else {
    chromosome.path <- FULL.PATH
  }
  
  # if not loaded yet
  if (is.null(genome[chromosome][[1]])) {

    chr.seq <- paste( scan(chromosome.path, "", skip = 1, quiet = TRUE), collapse = "")
    
    if (form == "vector") {
      
      chr.seq <- strsplit(chr.seq, split = "", fixed = TRUE)[[1]]
      
      if (vector.factor) {
        chr.seq <- as.factor(chr.seq)
      }
      
    }
    
    genome[chromosome] <<- list(chr.seq)
    
    cat(chromosome, "is loaded.\n")
    
  } else {
    cat(chromosome, "is already loaded.\n")
  }
}