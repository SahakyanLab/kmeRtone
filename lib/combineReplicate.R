combineReplicate <- function(list.genomic.coordinate) {
  # NOTE: havent been tested yet for big table
  
  # add column replicate_number
  cnt <- 1
  for (dt in list.genomic.coordinate) {
    dt[, replicate_number := cnt]
    cnt <- cnt + 1
  }
  
  dt <- rbindlist(list.genomic.coordinate)
  
  # detect replicate with the same coordinate and add column replicate
  #   - for more than one replicate, they will be pu together separated with space e.g. "1 2 3"
  dt[, replicate := paste(replicate_number, collapse = " "), by = .(chromosome, strand, start, end)]
  
  # remove column replicate_number
  dt[, replicate_number := NULL]
  
  # remove any duplicates
  dt <- unique(dt)
  
  return(dt)
}