readGenome <- function(chromosome, start, end, form="vector", case="upper",
         genome.path, prefix="", suffix=".fa.gz",
         FULL.PATH=NULL, size=F) {
  # This function...
  # - load on demand
  #    > only load a specific requested sequence
  # - it's core function, scan() become noticeably slow when requesting sequence
  #   at line >1,000,000 of fasta file
  #    > maybe Rcpp can speed things up - 
  # - will output character(0) if index out of range or NA if there is an empty last line
  #    > to give a friendly error require to calculate total number of lines,
  #      which at the moment is slow.
  #
  # Arguement
  # form       "vector" or "string"
  #            "string" support multiple start and end for vectorisation
  # case       "original", "upper", or "lower"
  
  # Dependencies: 
  
  if (is.null(FULL.PATH)) {
    # get chromosome path
    chromosome.path <- paste0(genome.path, prefix, chromosome, suffix) 
  } else {
    chromosome.path <- FULL.PATH
  }
  
  # get bases number per line
  bases.per.line <- nchar(scan(chromosome.path, "", skip = 1, nlines = 1, quiet = T))
  
  # calculate sequence length
  if (size == TRUE){
    
    # last line number
    last.line <- length(count.fields(chromosome.path, skip = 1))
    
    # number of bases at last line
    bases.at.last.line <- nchar(scan(chromosome.path, "", skip = last.line, nlines = 1, quiet = T))
    
    # length of sequence
    DNA.length <- (bases.per.line * (last.line-1)) + bases.at.last.line
    
    return(DNA.length)
  }
  
  # coordinate checking
  if (sum(start > end) > 0) {
    stop("End coordinate cannot be bigger than the start.")
  } else if (sum(start < 0 | end < 0) > 0) {
    stop("Out of range genomic coordinate.")
  }
  
  # get min start and max start number
  min.start <- min(start)
  max.end <- max(end)
  
  # which line does the start coordinate located
  line.start <- min.start %/% bases.per.line
  if ( min.start %% bases.per.line > 0 ) {
    line.start <- line.start + 1
  }
  
  # which line does the end coordinate located
  line.end <- max.end %/% bases.per.line
  if ( max.end %% bases.per.line > 0 ) {
    line.end <- line.end + 1
  }
  
  line.to.skip <- line.start - 1
  total.lines <- line.end - line.start + 1
  
  # fetch DNA segment which contains the requested DNA sequence
  # line to skip always plus one to skip a header
  dna.segment <- scan(chromosome.path, "", skip = line.to.skip+1, nlines = total.lines, quiet = T)
  
  # concatenate the lines
  dna.segment <- paste(dna.segment, collapse = "")
  
  # start and end index at dna segment
  idx.start <- start - (line.to.skip * bases.per.line)
  idx.end <- end - (line.to.skip * bases.per.line)
  
  # requested dna sequence
  dna.seq <- substring(dna.segment, idx.start, idx.end)
  
  # because there is a convention of upper and lowercase...
  if (case == "upper") {
    dna.seq <- toupper(dna.seq)
  } else if (case == "lower") {
    dna.seq <- tolower(dna.seq)
  }
  
  if (form == "vector") {
    dna.seq <- strsplit(dna.seq, "")
    if (length(dna.seq) == 1) dna.seq <- unlist(dna.seq)
  } 
  
  return(dna.seq)
}