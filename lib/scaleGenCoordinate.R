scaleGenCoordinate <- function(genomic.coordinate, genome, env=parent.frame(),
                               scale, side="both", na.omit=F) {
  # DESCRIPTION
  # Scale up or down the genomic coordinates. Any out of range numbers will be removed.
  #
  # genomic.coordinate   <string>    A data.table variable (to update by reference)
  # genome               <string>    A genome variable (not to copy in the local function environment)
  # env               <environment>  An environment where data.table and genome located.
  # scale                <integer>   Scale can be positive (scale up) or negative (scale down)
  # side                 <string>    Side to scale. Options are "upstream", "downstream", and "both"
  
  # Dependencies:
  #     Kmertone variables: genomic.coordinate, chromosome.names
  #     Packages          : data.table
  #     Functions         : readGenome, reverseComplement
  
  if (sum(class(genomic.coordinate) != "character") > 0 | sum(class(genome) != "character") > 0) {
    stop("Please input variable name in a string format.")
  }
  
  
  if (side == "both") {
    
    env[[genomic.coordinate]][ , start := start - scale]
    env[[genomic.coordinate]][ , end := end + scale]
    
  } else if (side == "upstream") {
    
    env[[genomic.coordinate]][strand == "+", start := start - scale] #$start <- genomic.coordinate$start - scale
    env[[genomic.coordinate]][strand == "-", end := end + scale] #$end <- genomic.coordinate$end + scale
    
  } else if (side == "downstream") {
    
    env[[genomic.coordinate]][strand == "+", end := end + scale] #$end <- genomic.coordinate$end + scale
    env[[genomic.coordinate]][strand == "-", start := start - scale] #$start <- genomic.coordinate$start - scale
    
  }
  
  # remove out of range genomic coordinate
  # less than genome chromosome start number i.e. < 1
  env[[genomic.coordinate]][start < 1, start := NA]
  env[[genomic.coordinate]][end < 1, end := NA]
  
  # # more than chromosome end number
  env[[genomic.coordinate]][start > attr(env[[genome]], "length")[chromosome], start:=NA]
  env[[genomic.coordinate]][end > attr(env[[genome]], "length")[chromosome], end:=NA]
  
  # cat("Below are the index numbers which have out of range coordinates:\n", 
  #     which(is.na(genomic.coordinate$start) | is.na(genomic.coordinate$end)))
  
  # remove the out of range number
  if (na.omit == TRUE) {
    env[[genomic.coordinate]] <- na.omit(env[[genomic.coordinate]], cols = c("start", "end"))
  }
}
