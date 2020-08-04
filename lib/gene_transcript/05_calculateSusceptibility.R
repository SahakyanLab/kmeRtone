calculateSusceptibility <- function(element.coordinate, genome, ref.kmers,
                                    cutoff, output) {
  # A gene can have multiple block of intron and CDS. Thus, susceptibility of
  # intron and CDS of a gene is calculated as overall susceptibility, not as
  # individual block. A gene can have two IGRs on it's right and left flanking
  # region. Susceptibility of IGR is calculated as individual block. Upstream
  # and downstream of a gene can overlap with neighbouring genes. Should I
  # remove them.

  transcript.susceptibility <-
    calculateTranscriptSusceptibility(element.coordinate,
                                      genome, ref.kmers,
                                      cutoff,
                                      output,
                                      verbose = TRUE,
                                      ncpu)

  assign("transcript.susceptibility", transcript.susceptibility,
         envir = parent.frame())

}
