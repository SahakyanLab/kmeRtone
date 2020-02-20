downloadMassGenome <- function(assembly.summary.path, folder.path) {
  # Automatically download COMPLETE GENOME of a SINGLE STRAIN for each species.
  # This faq website is referred: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq
  # Any special character ()[]'/ in asm_name is changed to _
  
  # Dependencies
  #   Packages: data.table, stringi
  
  if(file.exists(assembly.summary.path)){
    summary <- fread(assembly.summary.path)
  } else {
    
    # download assembly summary
    download.file(assembly.summary.path, destfile = paste0(folder.path, "assembly_summary.txt"))
    
    summary <- fread(paste0(folder.path, "assembly_summary.txt"), sep = "\t", quote = "")
  }
  
  dir.create(folder.path, showWarnings = F)

  # # assembly summary url
  # bacteria = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
  # viral.ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt"
  
  # only take latest version of complete genomes.
  summary <- summary[assembly_level == "Complete Genome" & version_status == "latest"]
  
  # take only the first strain of each organism
  idx <- summary[, !duplicated(stri_extract_first_regex(organism_name, ".*? [a-zA-Z]*"))]
  summary <- summary[idx]
  
  # save the filtered assembly summary for future reference
  fwrite(summary, paste0(folder.path, "filtered_assembly_summary.csv"))
  
  # change to https protocol
  summary[, ftp_path:=stri_replace_first_fixed(ftp_path, pattern = "ftp", replacement = "https")]
  
  # download the filtered genomes
  summary[, download.file(url = paste0(ftp_path, "/", `# assembly_accession`, "_",
                                       stri_replace_all_regex(asm_name, " |\\(|\\)|\\[|\\]|\\'|/", "_"),
                                       "_genomic.fna.gz"),
                          destfile = paste0(folder.path, `# assembly_accession`, ".fna.gz"),
                          quiet = T)
          ,
          by = 1:nrow(summary)]
}
