prepGenome <- function(genome.name, genome.path, genome.prefix,
                       genome.suffix, genome, env) {
  
  # Dependencies:
  #     Package: Biostrings
  #     
  
  if (class(genome) == "genome") {
    
    # add print function
    env$print.genome <- function(obj) print(attr(obj, "length"))
    env$genome <- genome
    
    cat("already loaded...")

  } else if (!is.null(genome.path)) {
    
    if(grepl("[RDSrds]", genome.path)){
      
      env$genome <- readRDS(genome.path)
      
      # add print function
      env$print.genome <- function(obj) print(attr(obj, "length"))
      
    } else {
      cat("\n -- Genome is saved as RDS file at:", paste0(genome.path, "/", genome.name, ".rds."))
      genome <- loadGenome(genome.path, genome.name, env, save.as.RDS = TRUE)
      cat("\n -- Loading genome...")
    }
    
  } else if (genome.name %in% c("hg19", "hg38")) {
    
    genome.path <- paste0("data/genome/human/", genome.name, "/", genome.name, ".rds")
    
    # Download and setup genome if not exist locally
    if(!file.exists(genome.path)){
      
      genome.folder <- paste0("data/genome/human/", genome.name, "/")
      suppressWarnings(dir.create(genome.folder, recursive = TRUE))
      
      # Boolean if fasta file exist in genome folder.
      fasta.not.exist <- length(list.files(genome.folder, pattern = "chr.+[Ff]a")) < 23
      
      # Download
      # The ftp storage hierachy is not consistent for hg19 and hg38. Some genome
      # is an old version. They are all over the place. I had to manually hunt down
      # for the url for the latest version.
      if(fasta.not.exist){
        
        # Final output of Download step is separate fasta files in genome folder
        
        cat("\n -- Genome is not available locally. Downloading...")
        
        if(genome.name == "hg38"){
          
          # This url link to a compressed folder "chroms" containing uncompressed fasta files.
          url.path <- "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chromFa.tar.gz"
          download.file(url.path, paste0(genome.folder, "hg38.chromFa.tar.gz"), quiet = TRUE)
          
          # Extract
          untar(paste0(genome.folder, "hg38.chromFa.tar.gz"), exdir = genome.folder)
          
          # Move file from genome/chroms folder to genome/
          chr.filenames <- list.files(paste0(genome.folder, "chroms/"))
          file.rename(paste0(genome.folder, "chroms/", chr.filenames), paste0(genome.folder, chr.filenames))
          
          # Delete chroms folder and hg38.chromFa.tar.gz
          unlink(paste0(genome.folder, "chroms"), recursive = TRUE)
          file.remove(paste0(genome.folder, "hg38.chromFa.tar.gz"))
          
        } else if(genome.name == "hg19"){
          
          genome.fasta.path <- paste0(genome.folder, genome.name, ".fa.gz")
          
          url.path <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
          download.file(url.path, genome.fasta.path, quiet = TRUE)
          
          splitFasta(genome.fasta.path, genome.folder, compression = 1)
          
          chr.filenames <- list.files(genome.folder, "chr[0-9A-Za-z_]+\\.(fa|fna|fasta)")
        }
      }
      
      cat("DONE!\n -- Building genome...")
      # Build genome
      env$genome <- loadGenome(genome.folder, genome.name, env, save.as.RDS = TRUE)
      cat("DONE!\n -- Loading genome...")
      
      # Final file organisation
      # Delete fasta files
      if("chr.filenames" %in% ls()) {
        file.remove(paste0(genome.folder, chr.filenames))
        suppressWarnings(file.remove(paste0(genome.folder, genome.name, ".fa.gz")))
      }
      # Write a text file containing record of current chromosomes in genome.RDS (chromosome name + length)
      chr.info <- data.table(chromosome = names(env$genome), length = attr(env$genome, "length"))
      fwrite(chr.info, paste0(genome.folder, "chr_list.tsv"), sep="\t")
    }
    
    env$genome <- readRDS(genome.path)
    
    # add print function
    env$print.genome <- function(obj) print(attr(obj, "length"))
    
  } else {
    
    stop(paste0("\nGenome ", genome.name, " is not available!"))
    
  }
}