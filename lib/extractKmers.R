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
  #     Package   : data.table
  #     Function  : reverseComplement
  #     Object    : <genome>
  
  if (!is.null(DNA.pattern)) {
    expansion.factor <- (k-nchar(DNA.pattern))/2
    pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern))
  }
  
  # remove region with lower than size k
  if (env1[[genomic.coordinate]][end - start + 1 < k, .N] > 0) {
    cat("--Removing regions with shorter than k length\n")
    env1[[genomic.coordinate]] <- env1[[genomic.coordinate]][end - start + 1 >= k] 
  }
  
  # To divide and process table by 100k rows - mayble help with memory efficiency.
  #    - Nothing changed for 13 million rows. Maybe above 13 million rows.
  #table.chunks <- gl(nrow(genomic.coordinate)/100000, 100000, nrow(genomic.coordinate))

  kmers <- env1[[genomic.coordinate]][, {
    
    # ----------------------------------------------------
    # DNA pattern specified
    if (!is.null(DNA.pattern)) {
      
      DNA.seq <- substring(env2[[genome]][[chromosome]], start, end)
      
      max.len <- max(nchar(DNA.seq))
      if (is.na(max.len)) stop("\nThere is NA in the table!")
      
      # sliding windows
      kmers <- substring(DNA.seq, 1:(max.len-k+1), k:max.len)
      
      if (strand == "-") kmers <- reverseComplement(kmers, form = "string")
      
      idx.case <- substring(kmers, pattern.pos[1], pattern.pos[length(pattern.pos)]) == DNA.pattern
      
      kmers <- kmers[idx.case]
      
      kmers <- table(kmers)
      
      # -------------------------------------------------
      # No DNA pattern specified
    } else if (is.null(DNA.pattern)) {

      DNA.seq <- substring(env2[[genome]][[chromosome]], start, end)
      
      max.len <- max(nchar(DNA.seq))
      if (is.na(max.len)) stop("There is NA in the table!")
      
      # sliding windows
      kmers <- substring(DNA.seq, 1:(max.len-k+1), k:max.len)
      
      if (strand == "-") kmers <- reverseComplement(kmers, form = "string")
      
      kmers <- table(kmers)

    }
    
    list(kmer = names(kmers), count = as.vector(kmers))
    
  }, by = .(chromosome, strand)][, .(kmer, count)]
  
  # aggregate the count
  kmers <- kmers[, list(count = sum(count)), by = kmer]
  
  # remove kmer with less than size k; this is a side effect of combining
  #   unequal size of region and vectorise it
  kmers <- kmers[nchar(kmer) == k]
  
  return(kmers)
}