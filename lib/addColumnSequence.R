addColumnSequence <- function(env) {
  # Because data.table is already multithreading, we can only see faster speed after 1M table row 
  # when using foreach loop
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, genome
  #     Packages          : data.table, (if ncpu > 1: foreach, doParallel)
  #     Function          : reverseComplement

  # sort
  setkey(env$genomic.coordinate, chromosome, start, end)
  
  # expand "*" to "+" and "-" if any
  if (nrow(env$genomic.coordinate[strand == "*"]) > 0) {
    message("WARNING! There is same genomic coordinate on plus and minus strands.")
    env$genomic.coordinate <- rbind(env$genomic.coordinate[strand %in% c("+", "-")],
                                    env$genomic.coordinate[strand == "*"][, strand := "+"],
                                    env$genomic.coordinate[strand == "*"][, strand := "-"])
    env$genomic.coordinate <- unique(env$genomic.coordinate[, .(chromosome, start, end, strand)])
                              
    gc()
  }
  
  if (ncpu == 1) {
    
    # divide table by 100000
    #table.chunks <- gl(nrow(env$genomic.coordinate)/100000, 100000, nrow(env$genomic.coordinate))
    
    # get sequence
    env$genomic.coordinate[, sequence := substring(env$genome[[chromosome]], start, end),
                       by = .(chromosome)]
    
    # divide table by 10000 for - strand
    #table.neg.chunks <- gl(nrow(env$genomic.coordinate[strand == "-"])/10000, 10000,
    #                       nrow(env$genomic.coordinate[strand == "-"]))
    
    # reverse complement for minus strand
    env$genomic.coordinate[strand == "-", sequence := reverseComplement(sequence, form = "string")]
    
  } else if (ncpu > 1) {
    
    # distribute rows to cpu
    chunks <- distributeChunk(nrow(env$genomic.coordinate), ncpu)
    
    genomic.coordinate.list <- lapply(1:ncpu, function(i) env$genomic.coordinate[ chunks$start[i]:chunks$end[i] ])
    
    # map the sequence
    toExport <- c("reverseComplement")
    seqs <- foreach(genomic.coordinate=genomic.coordinate.list, .combine = "c",
                    .noexport = ls()[!ls() %in% toExport], .packages = "data.table") %dopar% {
                      
      # get sequence
      genomic.coordinate[, sequence := substring(env$genome[[chromosome]], start, end),
                         by = .(chromosome)]
      
      # reverse complement for minus strand
      genomic.coordinate[strand == "-", sequence := reverseComplement(sequence, form = "string")]
      
      return(genomic.coordinate[, sequence])
    }
    
    env$genomic.coordinate[, sequence := seqs]
    
  }
}