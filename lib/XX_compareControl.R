compareControls <- function(kmers, genome.name, output=".", genome.kmers.path=NULL,
                            genome=NULL, genome.path=NULL, genome.prefix=NULL,
                            genome.suffix=NULL, genome=NULL, env=parent.frame()) {
  # Only support one size for k
  # 
  # kmers     <data.table>    data.table of kmers. It should contain 3 columns: kmer, case, control.
  #             <string>      path to kmers csv file.
  
  # Resolve kmers
  if(class(kmers)[1] == "character") kmers <- fread(kmers)
  
  k <- kmers[, unique(nchar(kmer))]
  
  # Resolve genome.kmers.path
  if(is.null(genome.kmers.path)) genome.kmers.path <- paste0("data/kmers/kmers_genome-", genome.name, ".csv")
  
  # Get genome kmers
  if(file.exists(genome.kmers.path)){
    
    genome.kmers <- fread(genome.kmers.path)
    
  } else{
    
    genome.kmers <- extractGenomeKmers(genome.name, k, genome, 
                                       genome.path, genome.prefix,
                                       genome.suffix, genome)
  }
  
  # Only take matching kmers
  genome.kmers <- genome.kmers[kmer %in% kmers$kmer]
  
  # Calculate proportion
  genome.kmers[, proportion := count/sum(count)]
  kmers[, proportion := control/sum(control)]
  
  # Sort kmer
  setkey(kmers, kmer)
  setkey(genome.kmers, kmer)
  
  # Plot
  pdf(paste0(output, "/control_comparison.pdf"))
  
  plot(genome.kmers$proportion, kmers$proportion, bty='l')
  abline(a = 0, b=1, col = "gray", lwd = 2, lty = 2)
  
  dev.off()
  
  # Calculate pearson correlation
  pearson.cor <- cor(genome.kmers$count, kmers$control, method = "pearson")
  
  # Perform chi square with p-value simulation
  chi2.test <- chisq.test(genome.kmers$count, kmers$control, simulate.p.value = TRUE)
  
  # Save correlation result
  writeLines(paste0("Pearson's correlation: ", pearson.cor), paste0(output, "/control_comparison.txt"))
  writeLines("\nChi square test with p-value simulation\n", paste0(output, "/control_comparison.txt"))
  writeLines(paste0("p-value: ", chi2.test$p.value), paste0(output, "/control_comparison.txt"))
  
  
}