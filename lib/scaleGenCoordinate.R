scaleGenCoordinate <- function(genomic.coordinate, chromosome.sizes, scale, side="both") {
  # DESCRIPTION
  # Scale up or down the genomic coordinates. Any out of range numbers will be removed.
  #
  # chromosomes.size       a named vector of chromosome size with chromosome names. e.g.
  #
  #                        > chromosome.sizes
  #                             chr1      chr2      chr3      chr4      chr5 ...
  #                        249250621 243199373 198022430 191154276 180915260 ...
  #
  # genomic.coordinate     a data.table; must have these columns: chromosome, start, end, strand
  # scale                  an integer - can be positive (scale up) or negative (scale down)
  
  # Dependencies: data.table
  
  if (side == "both") {
    
    genomic.coordinate[ , start := start - scale]
    genomic.coordinate[ , end := end + scale]
    
  } else if (side == "upstream") {
    
    genomic.coordinate[strand == "+", start := start - scale] #$start <- genomic.coordinate$start - scale
    genomic.coordinate[strand == "-", end := end + scale] #$end <- genomic.coordinate$end + scale
    
  } else if (side == "downstream") {
    
    genomic.coordinate[strand == "+", end := end + scale] #$end <- genomic.coordinate$end + scale
    genomic.coordinate[strand == "-", start := start - scale] #$start <- genomic.coordinate$start - scale
    
  }
  
  # remove out of range genomic coordinate
  # less than genome chromosome start number i.e. < 1
  genomic.coordinate[start < 1, start := NA]
  genomic.coordinate[end < 1, end := NA]
  
  # # more than chromosome end number
  
  genomic.coordinate[start > chromosome.sizes[chromosome], start:=NA, by = chromosome]
  genomic.coordinate[end > chromosome.sizes[chromosome], end:=NA, by = chromosome]
  
  # cat("Below are the index numbers which have out of range coordinates:\n", 
  #     which(is.na(genomic.coordinate$start) | is.na(genomic.coordinate$end)))
  
  # remove the out of range number 
  genomic.coordinate <- na.omit(genomic.coordinate, cols = c("start", "end"))
  
  return(genomic.coordinate)
}
