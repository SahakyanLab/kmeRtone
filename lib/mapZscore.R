mapZscore <- function(genome.name, kmers.ref, k, DNA.pattern) {
  # Method adopted from Claudia's thesis section 3.4 Cross-species Comparison.
  # 1. To calculate average z score of a given genome
  # 2. To plot fraction of DNA pattern out of whole genome vs. fraction of top 5% kmers
  #     a. To calculate fraction of top 5% kmers
  #     b. To calculate fraction of DNA pattern out of whole genome
  #
  # Output of this function is a vector containing average z score, fraction of top 5%,
  # and fraction of DNA pattern out of whole genome.
  #
  # genome.name   <string>      Genome name. List of genome available are hg19, hg38,
  #                             yokozuna1,
  # kmers.ref   <data.table>    Table of kmers of reference. One of the column should be
  #               <string>      named "z". It can a path to the table.
  # k             <numeric>     Kmer size
  # DNA.pattern   <string>      DNA pattern
  
  # load kmers reference i.e. with z score if path is provided
  if (class(kmers.ref)[1] == "character"){
    kmers.ref <- fread(kmers.ref)
  }
  
  genome.kmers.path <- paste0("data/kmers/", organism, "/kmers_k-", k, ".csv")
  genome.pattern.mer.path <- paste0("data/kmers/", organism, "/kmers_k-", nchar(DNA.pattern), ".csv")
  
  # load kmers to score
  genome.kmers <- fread(genome.kmers.path)
  
  # select idx of matched pattern of kmers.ref
  pattern.idx <- match(kmers.ref$kmer, genome.kmers$kmer)
  
  # assign z score
  genome.kmers[pattern.idx, z := kmers.ref$z]
  
  # calculate average z score
  average.z <- genome.kmers[pattern.idx, sum(count * z)/sum(count)]
  
  # identify cutoff of top 5% z score i.e. min z score of top 5%
  setorder(kmers.ref, -z)
  cutoff <- kmers.ref[round(0.05*.N),]$z
  
  # calculate fraction of top 5%
  top.fraction <- genome.kmers[pattern.idx, sum(count[z>=cutoff])/sum(count)]
  
  # pattern content in the genome (w/o flanking region)
  pattern.mers <- fread(genome.pattern.mer.path)
  pattern.fraction <- pattern.mers[, count[kmer == "CT"]/sum(count)]

  output <- c(average.z, top.fraction, pattern.fraction)
  names(output) <- c("average_z", "top_fraction", "pattern_fraction")
  
  return(output)
}