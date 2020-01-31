addColumnSequence <- function() {
  # Because data.table is already multithreading, we can only see faster speed after 1M table row 
  # when using foreach loop
  
  # Dependencies:
  #     Global variables: genomic.coordinate
  #     Packages        : data.table, foreach, doParallel
  #     Functions       : readGenome, reverseComplement
  
  # sort
  setkey(genomic.coordinate, chromosome, start, end)
  
  # expand "*" to "+" and "-" if any
  if (nrow(genomic.coordinate[strand == "*"]) > 0) {
    message("WARNING! There is same genomic coordinate on plus and minus strands.")
    message("         Update by reference won't work. You need to assign object.")
    genomic.coordinate <- rbind(genomic.coordinate[strand %in% c("+", "-")],
                genomic.coordinate[strand == "*"][, strand := "+"],
                genomic.coordinate[strand == "*"][, strand := "-"])
  }
  
  if (ncpu == 1) {
    
    # get sequence
    genomic.coordinate[, sequence := readGenome(chromosome, start, end, form = "string", genome.path = genome.path),
                       by = chromosome]
    
    # reverse complement for minus strand
    genomic.coordinate[strand == "-", sequence := reverseComplement(sequence, form = "string")]
    
  } else {
    
    # distribute rows to cpu
    chunks <- distributeChunk(nrow(genomic.coordinate), ncpu)
    
    genomic.coordinate.list <- lapply(1:ncpu, function(i) genomic.coordinate[ chunks$start[i]:chunks$end[i] ])
    
    # map the sequence
    
    seqs <- foreach(genomic.coordinate=genomic.coordinate.list, .combine = "c", .noexport = "genomic.coordinate",
                    .packages = "data.table") %dopar% {
                      
      # get sequence
      genomic.coordinate[, sequence := readGenome(chromosome, start, end, form = "string", genome.path = genome.path),
         by = chromosome]
      
      # reverse complement for minus strand
      genomic.coordinate[strand == "-", sequence := reverseComplement(sequence, form = "string")]
      
      return(genomic.coordinate[, sequence])
    }
    
    genomic.coordinate[, sequence := seqs]
    
  }
}