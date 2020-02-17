countReverseComplementKmers <- function(kmers, env=parent.frame()) {
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
  
  # Smartly count the reverse complement
  env[[kmers]][, rc := reverseComplement(kmer, form="string")]
  rc.count <- env[[kmers]][match(rc, kmer), count]
  
  # Because some reverse complement sequence are not exist in the plus strand,
  # match() results in NAs. Change NA to zero
  if (length(rc.count[is.na(rc.count)]) > 0) {
    rc.count[is.na(rc.count)] <- 0
  }
  
  # Update the count
  env[[kmers]][, count := count + rc.count]
  
  # Add the non-existing rc sequence
  if (length(rc.count[is.na(rc.count)]) > 0) {
    env[[kmers]] <- rbind(env[[kmers]], data.table(kmer = env[[kmers]][rc.count==0, rc], count = 1), fill=TRUE)
  }
  
  # Remove rc column
  env[[kmers]][, rc := NULL]
  
}