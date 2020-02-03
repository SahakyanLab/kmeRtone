scaleGenCoordinate <- function(genomic.coordinate, genome, scale, side="both") {
  # DESCRIPTION
  # Scale up or down the genomic coordinates. Any out of range numbers will be removed.
  #
  # genomic.coordinate   <string>    A data.table variable (to update by reference)
  # genome               <string>    A genome variable (not to copy in the local function environment)
  # scale                <integer>   Scale can be positive (scale up) or negative (scale down)
  # side                 <string>    Side to scale. Options are "upstream", "downstream", and "both"
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, chromosome.names
  #     Packages          : data.table
  #     Functions         : readGenome, reverseComplement
  
  genomic.coordinate <<- eval(parse(text = genomic.coordinate))
  genome <<- eval(parse(text = genome))
  
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
