addColumnSequence <- function(genomic.coordinate, genome, env=parent.frame()) {
  # Because data.table is already multithreading, we can only see faster speed after 1M table row 
  # when using foreach loop
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, genome
  #     Packages          : data.table, (if ncpu > 1: foreach, doParallel)
  #     Function          : reverseComplement

  start.time <- Sys.time()
  cat("Adding sequence to the table...")
  
  # sort
  setkey(env[[genomic.coordinate]], chromosome, start, end)
  
  # expand "*" to "+" and "-" if any
  if (nrow(env[[genomic.coordinate]][strand == "*"]) > 0) {
    message("WARNING! There is same genomic coordinate on plus and minus strands.")
    env[[genomic.coordinate]] <- rbind(env[[genomic.coordinate]][strand %in% c("+", "-")],
                                    env[[genomic.coordinate]][strand == "*"][, strand := "+"],
                                    env[[genomic.coordinate]][strand == "*"][, strand := "-"])
    env[[genomic.coordinate]] <- unique(env[[genomic.coordinate]][, .(chromosome, start, end, strand)])
                              
    gc()
  }
  
  # for shrinking uv data. comment out to be correct 
  # shrink start and end coordinates to point to their midpoint
  env[[genomic.coordinate]][, `:=`(
    
    start = {
      
      len <- end-start+1
      even <- len %% 2 == 0
      odd <- len %% 2 == 1
      
      # if even length
      start[even] = start + len/2 - 1
      # if odd length
      start[odd] = start + len%/%2
      start
    },
    
    end = {
      # if even length
      end[even] = end - len/2 + 1
      # if odd length
      end[odd] = end - len%/%2
      end
    }
  )]
  
  gc()
  
  # if (ncpu == 1) {
    
    # get sequence
    env[[genomic.coordinate]][, sequence := {

      dna.seq <- unlist(stri_sub_all(genome[[chromosome]], start, end))

      if (strand == "-") dna.seq <- reverseComplement(dna.seq, form = "string")
      dna.seq
      },by = .(chromosome, strand)]
    
  # } else if (ncpu > 1) {
  #   
  # 
  #   # distribute rows to cpu
  #   chunks <- distributeChunk(nrow(env[[genomic.coordinate]]), ncpu)
  #   
  #   genomic.coordinate.list <- lapply(1:ncpu, function(i) env[[genomic.coordinate]][ chunks$start[i]:chunks$end[i] ])
  #   
  #   # map the sequence
  #   toExport <- c("reverseComplement")
  #   seqs <- foreach(genomic.coordinate=genomic.coordinate.list, .combine = "c",
  #                   .noexport = ls()[!ls() %in% toExport], .packages = "data.table") %dopar% {
  #                     
  #   # divide table by 10mil
  #   #table.chunks <- gl(env[[genomic.coordinate]][, .N/5e+6], 5e+6, env[[genomic.coordinate]][, .N])
  #                     
  #   # get sequence
  #   env[[genomic.coordinate]][, sequence := {
  #     if (.N > 7e+5) {
  #       start <- split(start, ceiling(seq_along(start)/7e+5))
  #       end <- split(end, ceiling(seq_along(end)/7e+5))
  #     }
  #     dna.seq <- unlist(stri_sub_all(genome[[chromosome]], start, end))
  #     if (strand == "-") dna.seq <- reverseComplement(dna.seq, form = "string")
  #     
  #     dna.seq
  #      },by = .(chromosome, strand)]
  #     
  #     return(genomic.coordinate[, sequence])
  #   }
  #   
  #   env[[genomic.coordinate]][, sequence := seqs]
  # }
  gc()
  time.diff <- Sys.time() - start.time
  cat("DONE! ---", time.diff[1], attr(time.diff, "units"))
}