GCcontent <- function(env) {
  
  stop("NOT properly implemented yet.")
  
  width = c(100, 500, 1000, 5000, 10000, 50000, 100000)
  
  # load calculated value
  if (env$genome.name %in% c("hg19", "GRCh37")) {
    GC.percent <- fread("data/GRCh37_GC_percentage.csv")
  } else if (env$genome.name %in% c("hg38", "GRCh38")) {
    GC.percent <- fread("data/GRCh38_GC_percentage.csv")
  } else {
    # calculate
    GC.percent <- countN("genome", c("C", "G"), width, env$kmertone.env, ncpu)
    fwrite(GC.percent, "data/genome_GC_percentage.csv")
  }
  
  # calculate GC content of the case region
  
  
}