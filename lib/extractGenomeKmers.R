extractGenomeKmers <- function(genome.name, k, genome, 
                               genome.path=NULL, genome.prefix=NULL,
                               genome.suffix=NULL, env=parent.frame()) {
  # Extract kmers from genome
  
  # Dependencies:
  #     Packages: data.table, stringi
  #     Function: prepGenome, countReverseComplementKmers, reverseComplement
  
  # Prepare genome
  prepGenome(genome.name, genome.path, genome.prefix, genome.suffix, genome, env)
  
  # Only take base chr and exclude mitochondria. Some genome use roman numeral (IXV) e.g. yeast genome
  chr.names <- names(genome)[grep("chr[0-9XYIV]+$", names(genome))]
  
  # Build genome genomic.coordinate table
  genome.dt <- data.table(chromosome = chr.names,
                          start = 1,
                          end = attr(genome, "length")[chr.names],
                          strand = "*")
  
  # Extract kmers on + strand
  kmers <- extractKmers(genome.dt, genome, k)
  
  # Count kmers on - strand
  countReverseComplementKmers("kmers")
  
  # write to data folder
  suppressWarnings(dir.create("data/kmers", recursive = TRUE))
  fwrite(kmers, paste0("data/kmers/kmers_genome-", genome.name, ".csv"))
  
  return(kmers)
}



