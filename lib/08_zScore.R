zScore <- function(kmers, env=parent.frame()) {
  # Calculate Z score and update the table
  
  # kmers   <string>   A variable name pointing to kmers<data.table>. The
  
  env[[kmers]][, z := {
    
    # total case count
    total.case <- sum(case)
    
    # calculate proportion
    p.case <- case / total.case
    p.control <- control / sum(control)
    
    # From Claudia's code
    z <- sqrt( total.case ) * (p.case - p.control) / sqrt( p.control * (1 - p.control) )
    
    # I inferred based on the thesis equation 1
    z <- (case - control) / sqrt( control * (1 - p.control) )
    
    z
  }]
}