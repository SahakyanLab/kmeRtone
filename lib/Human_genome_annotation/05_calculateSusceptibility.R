calculateSusceptibility <- function(element.coordinate, genome, ref.kmers, score, cutoff, output){
  # A gene can have multiple block of intron and CDS. Thus, susceptibility of intron and CDS of a gene
  # is calculated as overall susceptibility, not as individual block. 
  # A gene can have two IGRs on it's right and left flanking region. Susceptibility of IGR is calculated
  # as individual block.
  # Upstream and downstream of a gene can overlap with neighbouring genes. Should I remove them.
  
  ref.kmers <- fread(ref.kmers)
  
  susceptibility.count <- calculateTranscriptSusceptibility(element.coordinate, genome, ref.kmers, score, cutoff,
                                                            output, verbose=TRUE, ncpu = ncpu)
  
  assign("susceptibility.count", susceptibility.count, envir = parent.frame())

}
