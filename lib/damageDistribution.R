damageDistribution <- function() {
  # venndiagram
  # any experimental bias
  
  # Dependencies:
  #     Global variables: genomic.coordinate, DNA.pattern
  #     Package         : data.table
  #     Function        : -
  
  options( scipen=999)
  
  # if there are replicates
  if(sum(grep("replicate", colnames(genomic.coordinate))) > 0) {
    
    replicate.names <- colnames(genomic.coordinate)[grep("replicate", colnames(genomic.coordinate))]
    
    total.dmg.reps <- sum(sapply(replicate.names, function(rep) {
      genomic.coordinate[eval(parse(text = rep)) > 0, .N]
    }))
    cat("\nTotal damage site:\t\t", total.dmg.reps)
  }
  
  # Check for duplicated record
  dup.cnt <- genomic.coordinate[duplicated(genomic.coordinate[ ,1:5]), .N]
  if (dup.cnt > 0) message(paste("WARNING! There are", dup.cnt, "duplicates!"))
  
  # Check damage pattern
  pattern.percent <- genomic.coordinate[, table(sequence) / .N * 100]
  
  if(sum(grep("replicate", colnames(genomic.coordinate))) > 0) {
    pattern.percent.by.replicate <- sapply(replicate.names, function(rep) {
      genomic.coordinate[ eval(parse(text = rep)) > 0 , table(sequence) / total.dmg.reps * 100]
    }) 
  }
  
  ## Print message
  cat("\nTotal unique damage site\t:", nrow(genomic.coordinate))
  cat("\nDamage pattern\t\t\t:", DNA.pattern, "\n")
  cat("\nThe percentage of unique damage site is: \n")
  print(pattern.percent)
  if (sum(!names(pattern.percent) %in% DNA.pattern) > 0) {message("WARNING! Unexpected damage pattern.")}
  else {cat("\n")}
  if(sum(grep("replicate", colnames(genomic.coordinate))) > 0) {
    cat("Percentage of damage occurence in each replicate (out of total damages)\n")
    print(pattern.percent.by.replicate)
  }
  
  # Check if damage happened at the same position of both strands
  idx.both.strands <- duplicated(genomic.coordinate[ ,1:3]) | duplicated(genomic.coordinate[ ,1:3], fromLast = T)
  dmg.both.strands.genomic.coordinate <- genomic.coordinate[idx.both.strands][order(chromosome, start, end)]
  
  ## Print message
  if (nrow(dmg.both.strands.genomic.coordinate) > 0) {
    message("\nWARNING! There are unique damages occur at the same position of both plus and minus strands.")
    print(dmg.both.strands.genomic.coordinate)
    
    # Count the percentage of this occurence
    dmg.both.strands.percent <- dmg.both.strands.genomic.coordinate[, table(sequence) / nrow(genomic.coordinate) * 100]
    by.replicate <- sapply(replicate.names, function(rep) {
      dmg.both.strands.genomic.coordinate[eval(parse(text = rep)) > 0, table(sequence)]
    })
    
    message(paste("\nTotal percentage (out of total unique damages) of such occurence is",
                  sum(dmg.both.strands.percent), "%"))
    cat("The percentage breakdown is as follow:\n")
    print( dmg.both.strands.percent )
    
    if(sum(grep("replicate", colnames(genomic.coordinate))) > 0) {
      # Count the percentage of this occurence in each replicate
      cat("\nThe percentage of replicates displaying this behaviour are as follow:\n")
      print(dmg.both.strands.cnt.by.rep)
      
      ## If there is a count mismatch between complementary sequences, that means the strands are mixed
      ## and matched to produce this observation (same genomic coordinates on both plus and minus strand)
      if ( sum(dmg.both.strands.cnt.by.rep[c("A"),] != dmg.both.strands.cnt.by.rep[c("T"),]) > 0 | 
           sum(dmg.both.strands.cnt.by.rep[c("C"),] != dmg.both.strands.cnt.by.rep[c("G"),]) > 0 ) {
        
        dmg.both.strands.cnt.by.rep <- sapply(replicate.names, function(rep) {
          idx <- duplicated(genomic.coordinate[ !is.na(eval(parse(text = rep))), 1:3]) | 
            duplicated(genomic.coordinate[ !is.na(eval(parse(text = rep))), 1:3],
                       fromLast = T)
          tab <- genomic.coordinate[ !is.na(eval(parse(text = rep))) ][idx, table(sequence) /  
                                                                         total.dmg.reps * 100]
          
          # fill missing sequence with zero so that sapply output a nice matrix instead of a list
          i <- !(c("A", "C", "G", "T") %in% names(tab))
          if (sum(i) > 0) {
            for (b in c("A", "C", "G", "T")[i]) {
              tab[[b]] <- 0
            }
            tab <- tab[order(names(tab))]
          }
          return(tab)
        })
        message("\nWARNING! There is a mix and match between strands of different replicates")
        cat("Here is the real percentage within each replicate where damage occurs at the same position of both strands.\n")
        print(dmg.both.strands.cnt.by.rep)
      }
    }
    
  }
}