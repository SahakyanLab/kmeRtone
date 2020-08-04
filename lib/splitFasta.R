splitFasta <- function(fasta, output="./", compression=9){
  
  if(length(fasta) == 1) big.fasta <- readLines(fasta)
  else big.fasta <- fasta
  
  header.idx <- grep("^>", big.fasta)
  
  special.characters <- "[]()<>:;,./\'\"*&^%$@?|#"
  special.characters <- strsplit(special.characters, "")[[1]]
  
  seq.line.starts <- header.idx + 1
  seq.line.ends <- c(header.idx[-1] - 1, length(big.fasta))
  
  for(i in 1:length(header.idx)){
    
    mini.fasta <- big.fasta[ header.idx[i]:seq.line.ends[i] ]
    
    mini.fasta.name <- big.fasta[header.idx][i]
    mini.fasta.name <- substring(mini.fasta.name, 2)
    
    if(sum(special.characters %in% strsplit(mini.fasta.name, "")[[1]]) > 0){
      mini.fasta.name <- paste0("seq", i)
    }
    
    mini.fasta.path <- paste0(output, "/", mini.fasta.name, ".fa.gz")
    
    gz <- gzfile(mini.fasta.path, "w", compression=compression)
    writeLines(mini.fasta, gz)
    close(gz)
    
  }
}