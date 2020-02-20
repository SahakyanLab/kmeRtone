getPQS <- function(genome) {
  # method adopted from "Machine learning model for sequence-driven
  # DNA G-quadruplex formation" under Method section named
  # "Putative quadriplex sequence (PQS) motifs." pg. 8
  # pseudo-regex: { G(3+) N(1-12) }(3+) G(3+)
  
  # Dependencies
  #     Packages: stringi, data.table

  pqs <- stri_locate_all_regex(genome, "(G{3,}.{1,12}?){3,}G{3,}")
  pqs <- lapply(seq_along(pqs), function(i) {
    
    data.table(chromosome = names(genome)[i], start = pqs[[i]][,"start"], end = pqs[[i]][, "end"],
               strand = "+")
    
  })
  pqs <- rbindlist(pqs)
  
  # on minus strand
  pqs.rc <- stri_locate_all_regex(genome, "(C{3,}.{1,12}?){3,}C{3,}")
  pqs.rc <- lapply(seq_along(pqs.rc), function(i) {
    
    data.table(chromosome = names(genome)[i], start = pqs.rc[[i]][,"start"], end = pqs.rc[[i]][, "end"],
               strand = "-")
    
  })
  pqs.rc <- rbindlist(pqs.rc)
  
  # combine plus and minus strand
  pqs.dt <- na.omit(rbind(pqs, pqs.rc))
  
  return(pqs.dt)
}
