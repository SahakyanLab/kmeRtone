extractKmers <- function(genomic.coordinate, genome, k, DNA.pattern,
                         env=parent.frame()) {
  # Extract kmers from genomic coordinate table
  #
  # genomic.coordinate   <string>     A variable name pointing to genomic coordinate
  #                                   <data.table>. Strand can be +- or *. * strand
  #                                   is treated as + strand only.
  # genome               <string>     A <genome> class object.
  # env                <environment>  An enviroment object for genomic coordinate.
  # k                    <numeric>    Size of kmer. <integer> is accepted as well.
  # DNA.pattern          <string>     DNA pattern at the center of the kmers.
  
  # Dependencies
  #     Package   : data.table, stringi
  #     Function  : reverseComplement
  #     Object    : <genome>
  
  # calculate expandion factor and pattern position
  if (!is.null(DNA.pattern)) {
    expansion.factor <- (k-nchar(DNA.pattern))/2
    pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern))
  }
  
  # all possible kmers - these are used later for fast binary matching %in%
  # This is also used for all kmers w/o pattern to remove base N
  possible.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
  if (!is.null(DNA.pattern)) {
    pattern.idx <- possible.kmers[, do.call(paste0,.SD) %in% DNA.pattern,
                                  .SDcols = pattern.pos]
    possible.kmers <- possible.kmers[pattern.idx]
  }
  possible.kmers <- possible.kmers[, do.call(paste0,.SD)]
  
  # get kmers
  kmers <- env[[genomic.coordinate]][!is.na(end-start), {
    
    mini.table.chunk <- distributeChunk2(end-start+1-k+1, 1e+7, "label")
      
      kmers <- .SD[, {
        
        if (sum(unique(end-start+1) != k) > 0) {
          start = lapply(1:.N, function(i) start[i]:(end[i]-k+1))
        }
        
        DNA.seq <- unlist(stri_sub_all(genome[[chromosome]], from = start, length = k))
        
        if (strand == "-") DNA.seq <- reverseComplement(DNA.seq, form = "string")
        
        DNA.seq <- DNA.seq[DNA.seq %in% possible.kmers]
        
        DNA.seq <- table(DNA.seq)
        
        list(kmer = names(DNA.seq), count = as.vector(DNA.seq))
        
      }, by = mini.table.chunk][, .(kmer, count)]
      
      kmers <- kmers[, list(count = sum(count)), by = kmer]
    
    kmers
  }, by = .(chromosome, strand)][, .(kmer, count)]
  
  # aggregate the count
  kmers <- kmers[, list(count = sum(count)), by = kmer]
  
  return(kmers)
}