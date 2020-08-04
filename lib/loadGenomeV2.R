loadGenome <- function(genome, chromosome, genome.path, genome.prefix, genome.suffix,
                       full.path=NULL, env=parent.frame(), form="string",
                       letter.case="upper") {
  
  # genome          <string>     A variable name pointing to a list object containing chromosome sequence.
  # chromosome      <string>     A chromosome name.
  # genome.path     <string>     A path to a folder containing chromosome fasta files.
  # genome.prefix   <string>     A prefix name of the fasta filename before the name of the chromosome.
  # genome.suffix   <string>     A suffix name of the fasta filename after the name of the chromosome.
  # full.path       <string>     One can provide a full path to the chromosome fasta file.
  # env           <environment>  An environment object where the genome object exists.
  # form            <string>     Output the sequence in either in a "vector" or "string" format.
  # letter.case     <string>     Output the sequence in either in an "upper" case or a "lower" case.
  
  start.time <- Sys.time()
  
  if (class(genome)[1] != "character") {
    stop("Please input genome variable name in a string format.")
  }
  
  if (!genome %in% ls(envir = env)) {
    
    env[[genome]] <- list()
    
    class(env[[genome]]) <- "genome"
    
    # add print function
    env$print.genome <- function(obj) print(attr(obj, "length"))
  }
  
  if (is.null(full.path)) {
    # get chromosome path
    chromosome.path <- paste0(genome.path, genome.prefix, chromosome, genome.suffix)
    
  } else {
    chromosome.path <- full.path
  }
  
  # if not loaded yet
  if (is.null(env[[genome]][chromosome][[1]])) {
    
    chr.seq <- paste( scan(chromosome.path, "", skip = 1, quiet = TRUE), collapse = "")
    
    if (form == "vector") {
      
      chr.seq <- strsplit(chr.seq, split = "", fixed = TRUE)[[1]]
      
      if (vector.factor) {
        chr.seq <- as.factor(chr.seq)
      }
      
    }
    
    env[[genome]][chromosome] <- list(chr.seq)
    
    # add length as attribute
    attr(env[[genome]], "length") <- c(attr(env[[genome]], "length"), chr = nchar(chr.seq))
    names(attr(env[[genome]], "length"))[names(attr(env[[genome]], "length")) == "chr"] <- chromosome
    
    time.diff <- Sys.time() - start.time
    cat(chromosome, "is loaded. ---", time.diff[1], attr(time.diff, "units"), "\n")
    
  } else {
    cat(chromosome, "is already loaded.\n")
  }
  
}
