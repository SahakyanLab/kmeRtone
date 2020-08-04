extractKmers <- function(genomic.coordinate, genome, k, DNA.pattern,
                         env1=parent.frame(), env2=env1) {
  # Extract kmers from genomic coordinate table
  #
  # genomic.coordinate   <string>     A variable name pointing to genomic coordinate
  #                                   <data.table>. Strand can be +- or *. * strand
  #                                   is treated as + strand only.
  # genome               <string>     A variable name for <genome> class object.
  # env1               <environment>  An enviroment object for genomic coordinate.
  # env2               <environment>  An enviroment object for genome.
  # k                    <numeric>    Size of kmer. <integer> is accepted as well.
  # DNA.pattern          <string>     DNA pattern at the center of the kmers.
  
  # Dependencies
  #     Package   : data.table, stringi
  #     Function  : reverseComplement
  #     Object    : <genome>
  
  if (!is.null(DNA.pattern)) {
    expansion.factor <- (k-nchar(DNA.pattern))/2
    pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern))
  }
  
  # calculate expandion factor and pattern position
  expansion.factor <- (k-nchar(DNA.pattern))/2
  pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern))
  
  # initiate all possible kmers
  # initiate count; should be rep(NA, all.possible.kmer)
  possible.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
  pattern.idx <- possible.kmers[, do.call(paste0,.SD) %in% DNA.pattern,
                       .SDcols = pattern.pos]
  possible.kmers <- possible.kmers[pattern.idx]
  possible.kmers <- possible.kmers[, do.call(paste0,.SD)]
  
  # remove region with lower than size k
  if (env1[[genomic.coordinate]][end - start + 1 < k, .N] > 0) {
    cat("--Removing regions with shorter than length k.\n")
    env1[[genomic.coordinate]] <- env1[[genomic.coordinate]][end - start + 1 >= k]
    gc()
  }
  
  # To divide and process table by 100k rows - to help with memory efficiency.
  #table.chunks <- gl(env1[[genomic.coordinate]][!is.na(end-start), .N]/100000, 100000,
  #                   env1[[genomic.coordinate]][!is.na(end-start), .N])
  #table.chunks <- 1:nrow(env1[[genomic.coordinate]][!is.na(end-start)])

  kmers <- env1[["genomic.coordinate"]][!is.na(end-start), {
    
    # ----------------------------------------------------
    # DNA pattern specified
    if (!is.null(DNA.pattern)) {
      
      
      count <- rep(0, length(possible.kmers))
      i <- 1
      while(i <= .N) {
        
        DNA.seq <- substring(env2[["genome"]][[chromosome]], start[i], end[i])
        if (strand == "-") DNA.seq <- reverseComplement(DNA.seq, form = "string")
        
        count <- count + stri_count_fixed(DNA.seq, possible.kmers, overlap=TRUE)
        
        i<-i+1
      }
      
      count
      
      #alternative
      # kmers <- lapply(1:.N, function(i) {
      #   
      #   DNA.seq <- substring(env2[["genome"]][[chromosome]], start[i], end[i])
      #   if (strand == "-") DNA.seq <- reverseComplement(DNA.seq, form = "string")
      #   
      #   DNA.seq <- strsplit(DNA.seq, "")[[1]]
      #   
      #   count <- rep(0, length(DNA.seq)-k+1)
      #   
      #   j<-1
      #   while (j <= length(DNA.seq)-k+1) {
      #     
      #     dna.seq <- DNA.seq[j:(j+k-1)]
      #     
      #     if (paste(dna.seq[pattern.pos], collapse="") %in% DNA.pattern) {
      #       count[j] <- count[j] + 1
      #       names(count)[j] <- paste(dna.seq, collapse = "")
      #     }
      #     
      #     j<-j+1
      #   }
      #   
      #   count <- count[count > 0]
      #   return(count)
      # })
      # kmers <- unlist(kmers)
      # kmers <- tapply(kmers, names(kmers), sum)
      # 
      # kmers
      
      # -------------------------------------------------
      # No DNA pattern specified
    } else if (is.null(DNA.pattern)) {

      DNA.seq <- substring(env2[[genome]][[chromosome]], start, end)
      
      max.len <- max(nchar(DNA.seq))
      if (is.na(max.len)) stop("\nThere is NA in the table!")
      
      # sliding windows
      kmers <- stri_sub_all(DNA.seq, 1:(max.len-k+1), k:max.len)
      kmers <- unlist(kmers)
      kmers <- kmers[kmers != ""]
      
      if (strand == "-") kmers <- reverseComplement(kmers, form = "string")
      
      kmers <- table(kmers)

    }
    
    list(kmer = possible.kmers, count = count)
    
  }, by = .(chromosome, strand)][, .(kmer, count)]
  
  # aggregate the count
  kmers <- kmers[, list(count = sum(count)), by = kmer]
  gc()
  
  return(kmers)
}