countReverseComplementKmers <- function(kmers, env=parent.frame(), update.count=TRUE) {
  # Smartly count the reverse complement sequence from it's contrast strand
  #
  # kmers     <string>      A variable name pointing to data.table object. The data.table
  #                         must contains two named-columns: kmer and count
  # env     <environment>   An environment object where the data.table exist.
  
  # Dependency
  #     Package   : data.table, stringi
  #     Function  : reverseComplement
  
  if (class(kmers)[1] != "character") {
    stop("Please input a variable name of kmers table in string format.")
  }
  
  # Find reverse complement
  env[[kmers]][, rc.seq := reverseComplement(kmer, form="string")]
  
  # Some reverse complement sequence may not exist in the plus strand, so add the new sequence if any
  if(env[[kmers]][!rc.seq %in% kmer, .N > 0]){
    env[[kmers]] <- rbind(env[[kmers]], env[[kmers]][!rc.seq %in% kmer, .(kmer = rc.seq, count = 0, rc.seq = kmer)])
  }

  # Count complementary sequence
  rc.count <- env[[kmers]][match(rc.seq, kmer), count]
  
  # Because some reverse complement sequence are not exist in the plus strand,
  # match() results in NAs. Change NA to zero
  if (length(rc.count[is.na(rc.count)]) > 0) {
    stop("Something is wrong!")
    rc.count[is.na(rc.count)] <- 0
  }
  
  if(update.count){
    
    # Update the count
    env[[kmers]][, count := count + rc.count]
    
  } else {
    
    # Add column rc.count
    env[[kmers]][, rc.count := rc.count]
    
  }
  
  # Remove rc.seq column
  env[[kmers]][, rc.seq := NULL]
  
}