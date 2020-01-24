reverseComplement <- function(DNA.sequence) {
    # DNA.sequence must be in vector format e.g. c("A", "C", "G", "T", ...)
    
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
    
    return(DNA.sequence)
}