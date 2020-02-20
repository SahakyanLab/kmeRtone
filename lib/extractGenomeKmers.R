extractGenomeKmers <- function(genome, single.fasta.path=NULL, k) {
  # Extract kmers from genome
  # Notes: oligonucleotideFrequency is not stable at k > 10 for hg19
  #        and k > 11 for tardigrade genome
  
  # Dependencies:
  #     Packages: data.table, stringi, Biostrings
  #     Function: countReverseComplementKmers, reverseComplement
  
  if(k < 11){
    
    # load genome fasta
    if(!is.null(single.fasta.path)) {
      genome <- readDNAStringSet(fasta.path, use.names = F)
      
      # convert to DNAStringSet
    } else if(class(genome) %in% c("genome", "list")) {
      genome <- DNAStringSet(genome)
    }
    
    # extract kmers
    kmers <- oligonucleotideFrequency(genome, k)
    
    # add every column for every scaffold/ chromosome
    kmers <- colSums(kmers)
    
    # change to data.table
    kmers <- data.table(kmer = names(kmers), count = kmers)
    
  } else if (k > 10){
    
    cat("k is bigger than 10. It will take a long time.\n")
    
    # initiate all possible kmers
    kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
    kmers <- kmers[, .(kmer = do.call(paste0,.SD))]
    
    # count the kmers
    kmers[, count := {
      
      count <- sapply(genome, stri_count_fixed, pattern=kmer, overlap=TRUE)
      count <- rowSums(count)
      
      count
    }]
    
  }
  
  # count the opposite strand kmers
  countReverseComplementKmers("kmers")
  
  return(kmers)
}



