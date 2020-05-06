processGenomeAnnotationTable <- function(genome.annotation, genome, igr.rel.pos, igr.min.length,
                                         gene.upstream.length, gene.downstream.length){
  # AIM: To generate a genomic coordinate table for transcription sites. These elements are included
  # in the coordinate table: exon, intron, 5'-UTR, 3'-UTR, 2500 nt upstream, and 1000 nt downstream.
  #
  
  # Dependency 
  #    Package    : data.table
  #    Functions  : reverseComplement, processUCSCannoTable
  
  # If other than RefSeq annotation table, stop because not supported yet.
  if(!"cdsStartStat" %in% names(genome.annotation)) {
    stop("Sorry, this annotation table is not supported yet. Please contact developer for more support. Thank you.")
  }
  
  # Only take base chromosome and not mitochondria
  genome.annotation <- genome.annotation[chrom %in% c(paste0("chr", c(1:22, "X", "Y")))]
  
  # Create gene regions
  transcript.regions <- genome.annotation[, .(chromosome=chrom, start = txStart+1, end=txEnd, strand)]
  mergeGenCoordinate("transcript.regions")

  # Create intergenic region coordinates
  intergenic.coordinate <- transcript.regions[, .(chromosome,
                                                  start = c(start-igr.rel.pos[2], end+igr.rel.pos[1]),
                                                  end = c(start-igr.rel.pos[1], end+igr.rel.pos[2]),
                                                  strand = strand)]
  trimGenCoordinates("intergenic.coordinate", genome, remove=TRUE)
  mergeGenCoordinate("intergenic.coordinate")
  removeCaseZone("intergenic.coordinate", "transcript.regions", "sensitive")
  intergenic.coordinate <- intergenic.coordinate[end-start+1 >= igr.min.length]
  intergenic.coordinate[, c("name", if("name2" %in% colnames(genome.annotation)) "name2") := paste0("IGR_", 1:.N)]
  intergenic.coordinate[, element := "IGR"]
  
  rm(transcript.regions)
  
  # Only take complete coding sequence status
  genome.annotation <- genome.annotation[cdsStartStat == "cmpl" & cdsEndStat == "cmpl"]
  
  # Only take protein-coding transcript tagged as NM_ (NR_ is non-protein-coding transcript)
  genome.annotation <- genome.annotation[grepl("^NM_", name)]
  
  # For duplicate name2 (multi-splice transcript?), take the longest one.
  genome.annotation <- genome.annotation[order(txEnd-txStart, decreasing = TRUE)][!duplicated(name2)]
  
  # Extract elements of interest: CDS, intron, 5'-UTR, 3'-UTR, 2500 upstream, and 1000 nt downstream.
  element.coordinate <- processUCSCannoTable(genome.annotation, type = "refseq",
                                             element = c("cds", "intron", "utr"),
                                             gene.upstream.length, gene.downstream.length, genome)
  
  # Only take canonical gene transcript i.e. must have both 5'-UTR and 3'-UTR
  #  1) 5'-UTR, CDS, intron, 3'-UTR
  #  2) 5'-UTR, CDS, 3'-UTR
  setkey(element.coordinate, name)
  element.coordinate <- element.coordinate[, if(sum(unique(element) %in% c('UTR5', "UTR3")) == 2) .SD, by = name]
  
  assign("element.coordinate", rbind(element.coordinate, intergenic.coordinate), envir = parent.frame())
  #assign("intergenic.coordinate", intergenic.coordinate, envir = parent.frame())
  
}