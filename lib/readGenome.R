readGenome <- function(genome.path, chromosome, start, end, suffix=".fa.gz",
                       form="vector", case="upper") {
  # THE COMPARISON
  # GEN_loadGenome cannot accepts .gz file which load way faster than the
  # uncompressed counterpart.
  # Reason: GEN_loadGenome depends on UTIL_readLinesFast which uses readChar().
  #         The size of compressed file will underestimate the number of
  #         character for the readChar function.
  # This readGenome() function uses a bit slower function (scan vs. readChar)
  # to read the file but the compressed file makes it ~3x faster than GEN_loadGenome.
  # I expect by using Rcpp, this function can perform even much better.
  # 
  # This function...
  # - load on demand
  #    > only load a specific requested sequence
  # - it's core function, scan() become noticeably slow when requesting sequence
  #   at line >1e+6 of fasta file
  #    > maybe Rcpp can speed things up
  # - will output character(0) if index out of range
  #    > to give a friendly error require to calculate total number of lines,
  #      which at the moment is slow.
  #
  # Arguement
  # form       "vector" or "string"
  #             "string" support multiple start and end for vectorisation
  # case       "original", "upper", or "lower"
  
  # Dependency: stringi - for speed
  require(stringi)
  
  # get chromosome path
  chromosome.path <- paste0(genome.path, chromosome, suffix)

  # get bases number per line
  bases.per.line <- nchar(scan(chromosome.path, "", skip = 1, nlines = 1, quiet = T))
  
  # # calculate sequence length
  # ## last line number
  # last.line <- countLines(chromosome.path)[1] # noticably SLOW!
  # 
  # ## number of bases at last line
  # bases.at.last.line <- nchar(scan(chromosome.path, "", skip = last.line-1, nlines = 1, quiet = T))
  # 
  # ## length of sequence
  # DNA.length <- (bases.per.line * last.line-1) + bases.at.last.line
  
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
  
  # fetch DNA segment which contains the resuested DNA sequence
  # line to skip always plus one to skip a header
  dna.segment <- scan(chromosome.path, "", skip = line.to.skip+1, nlines = total.lines, quiet = T)
  
  if (form == "vector") {
    
    dna.segment <- unlist( strsplit(dna.segment, "") )
    
  } else if (form == "string") {

    dna.segment <- paste(dna.segment, collapse = "")
    
  }
  
  # start and end index at dna segment
  idx.start <- start - (line.to.skip * bases.per.line)
  idx.end <- end - (line.to.skip * bases.per.line)
  
  # requested dna sequence
  
  if (form == "vector") {
    
    dna.seq <- dna.segment[ idx.start:idx.end ]
    
  } else if (form == "string") {
    
    dna.seq <- stri_sub(dna.segment, idx.start, idx.end)

  }
  
  # because there is a convention of upper and lowercase...
  if (case == "upper") {
    dna.seq <- stri_trans_toupper(dna.seq)
  } else if (case == "lower") {
    dna.seq <- stri_trans_tolower(dna.seq)
  }
  
  return(dna.seq)
}