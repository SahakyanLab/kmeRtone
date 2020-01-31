reverseComplement <- function(DNA.sequence, form="vector") {
  # DNA.sequence can be a string or in a vector
  # form    "string" or "vector"
  #         "string" support vectorisation. e.g. c("AA", "TAATA")
  
  # Dependency: stringi (for fast vectorisation)
  
  if (form == "string") {

    suppressPackageStartupMessages( require(stringi) )
    
    # complement
    # equivalent to base function, chartr but a bit faster
    DNA.sequence <- stri_trans_char(DNA.sequence, "ACGT", "TGCA")
    
    # reverse
    DNA.sequence <- stri_reverse(DNA.sequence)
    
  } else if (form == "vector") {
    
    # locate each nucleotide base in the sequence
    idx.A <- DNA.sequence == "A"
    idx.C <- DNA.sequence == "C"
    idx.G <- DNA.sequence == "G"
    idx.T <- DNA.sequence == "T"
    
    # complement the sequence
    DNA.sequence[idx.A] <- "T"
    DNA.sequence[idx.C] <- "G"
    DNA.sequence[idx.G] <- "C"
    DNA.sequence[idx.T] <- "A"
    
    # reverse the sequence
    if (is.null(dim(DNA.sequence))) {
      DNA.sequence <- rev(DNA.sequence)
    } else {
      # if input is a matrix where each column is one DNA
      DNA.sequence <- apply(DNA.sequence, 2, rev)
    }
    
  }
  
  
  
  return(DNA.sequence)
}