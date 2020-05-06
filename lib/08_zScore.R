zScore <- function(kmers, env=parent.frame()) {
  # Calculate Z score and update the table
  # Claudis's thesis equation 1
  
  # kmers   <string>   A variable name pointing to kmers<data.table>. The
  
  env[[kmers]][, z := {
    
    # total case count (n)
    total.case <- sum(case)
    
    # proportion control (p)
    p.control <- control / sum(control)
    
    # predicted case distribution (np)
    case.predict <- total.case * p.control
    
    z <- (case - case.predict) / sqrt( case.predict * (1 - p.control) )
    
    z
  }]
}