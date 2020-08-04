extractKmers <- function(genomic.coordinate, genome, k, DNA.pattern=NULL) {
  # Extract kmers from genomic coordinate table
  #
  # genomic.coordinate   <data.table> Genomic coordinate
  #                                   Strand can be +, - or *. * strand
  #                                   is treated as + strand only.
  # genome               <string>     A <genome> class object.
  # k                    <numeric>    Size of kmer. <integer> is accepted as well.
  # DNA.pattern          <string>     DNA pattern at the center of the kmers.
  
  # Dependencies
  #     Package   : data.table, stringi
  #     Function  : reverseComplement, distributeChunk2
  #     Object    : <genome>
  
  # Input checking
  if(class(genome)[1] != "genome"){
    stop("Please input <genome> object for the genome.")
  } else if (class(genomic.coordinate)[1] == "character"){
    genomic.coordinate <- fread(genomic.coordinate, showProgress = FALSE)
  } else if (class(genomic.coordinate)[1] != "data.table"){
    stop("Please input genomic coordinate as <data.table> object.")
  }
  
  # Check genomic.coordinate - return empty data.table
  invalid.coordinate <- genomic.coordinate[!is.na(end-start), sum(end-start+1 >= k)] < 1
  if(invalid.coordinate) return(data.table(kmer=character(), count=numeric()))
  
  # all possible kmers - these are used later for fast binary matching %in%
  # This is also used for all kmers w/o pattern to remove base N
  
  if (!is.null(DNA.pattern)) {
    
    # calculate expansion factor and pattern position
    expansion.factor <- (k-nchar(DNA.pattern[1]))/2
    #pattern.pos <- seq(expansion.factor + 1, expansion.factor + nchar(DNA.pattern[1]))
    
    # include complementary sequence as well
    rc.DNA.pattern <- reverseComplement(DNA.pattern, form = "string")
    
    expansion.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), expansion.factor))
    expansion.kmers <- expansion.kmers[, do.call(paste0, .SD)]
    expansion.kmers.combn <- do.call(CJ, rep(list(expansion.kmers), 2))
    
    possible.kmers <- lapply(c(DNA.pattern, rc.DNA.pattern), function(dna.pattern) {
      
      possible.kmers <- expansion.kmers.combn[, .(V1, dna.pattern, V2)]

      return(possible.kmers)
    })
    
    possible.kmers <- rbindlist(possible.kmers)

  } else if(is.null(DNA.pattern)) {
    # k is limited to 15 because vector size is limited to .Machine$integer.max
    possible.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
  }
  
  possible.kmers <- possible.kmers[, do.call(paste0,.SD)]
  
  # get kmers
  kmers <- genomic.coordinate[ (!is.na(end-start)) & ((end-start+1) >= k), {
    
    mini.table.chunk <- distributeChunk2(end-start+1-k+1, 1e+7, "label")
    
    kmers <- .SD[, {
      
      if (sum(unique(end-start+1) != k) > 0) {
        start = lapply(1:.N, function(i) start[i]:(end[i]-k+1))
      }

      kmers <- unlist(stri_sub_all(genome[[chromosome]], from = start, length = k))
      
      if (strand == "-") kmers <- reverseComplement(kmers, form = "string")
      
      kmers <- kmers[kmers %in% possible.kmers]
      
      kmers <- table(kmers)
      
      counts <- as.vector(kmers)
      kmers <- names(kmers)
      
      # In an event where there is no kmer
      if(is.null(kmers)) kmers <- count <- NA
      
      list(kmer = kmers, count = counts)
      
    }, by = mini.table.chunk][, .(kmer, count)]
    
    kmers <- kmers[, list(count = sum(count)), by = kmer]
    
    kmers
  }, by = .(chromosome, strand)]
  
  kmers[, c("chromosome", "strand") := NULL]
  
  #print(kmers)
  
  # Omit NA (non-existence kmers i.e. bases other than ACTG)
  kmers <- na.omit(kmers)
  
  # aggregate the count
  kmers <- kmers[, list(count = sum(count)), by = kmer]
  
  return(kmers)
}