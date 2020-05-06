caseDistribution <- function(genomic.coordinate, DNA.pattern, output, env=parent.frame()) {
  # venndiagram
  # any experimental bias
  
  # Dependencies:
  #     Global variables: genomic.coordinate, DNA.pattern
  #     Package         : data.table
  #     Function        : -
  
  options( scipen=999)
  
  # replicate exist?
  replicate.exist <- sum(grep("replicate", colnames(env[[genomic.coordinate]]))) > 0
  
  # if there are replicates
  if(replicate.exist) {
    
    replicate.names <- colnames(env[[genomic.coordinate]])[grep("replicate", colnames(env[[genomic.coordinate]]))]
    
    total.case.reps <- sum(sapply(replicate.names, function(rep) {
      env[[genomic.coordinate]][eval(parse(text = rep)) > 0, .N]
    }))
    cat("Total case site\t\t\t:", total.case.reps, "\n")
  }
  
  # Check for duplicated record
  dup.cnt <- env[[genomic.coordinate]][duplicated(env[[genomic.coordinate]][ ,1:5]), .N]
  if (dup.cnt > 0) message(paste("WARNING! There are", dup.cnt, "duplicates!"))
  
  # Check case pattern
  if (!is.null(DNA.pattern)) {
    pattern.percent <- env[[genomic.coordinate]][, table(sequence) / .N * 100]
  }
  
  if(replicate.exist) {
    pattern.percent.by.replicate <- sapply(replicate.names, function(rep) {
      env[[genomic.coordinate]][ eval(parse(text = rep)) > 0 , table(sequence) / total.case.reps * 100]
    }) 
  }
  
  total.union <- nrow(env[[genomic.coordinate]])
  ## Print message
  cat("Total unique case site\t\t:", total.union, "\n")
  if (!is.null(DNA.pattern)) cat("\nDNA pattern\t\t\t:", DNA.pattern, "\n")
  
  cat("\nThe percentage of unique case site is: \n")
  print(pattern.percent)
  if (sum(!names(pattern.percent) %in% DNA.pattern) > 0) {message("WARNING! Unexpected DNA pattern.")}
  else {cat("\n")}
  
  if(replicate.exist) {
    cat("Percentage of case occurence in each replicate (out of total case)\n")
    print(pattern.percent.by.replicate)
    cat("\n")
  }
  
  # Check if case happened at the same position of both strands
  idx.both.strands <- duplicated(env[[genomic.coordinate]][ ,1:3]) | duplicated(env[[genomic.coordinate]][ ,1:3], fromLast = T)
  case.both.strands.genomic.coordinate <- env[[genomic.coordinate]][idx.both.strands][order(chromosome, start, end)]
  
  ## Print message
  if (nrow(case.both.strands.genomic.coordinate) > 0) {
    message("WARNING! There are unique cases occur at the same position of both plus and minus strands.")
    print(case.both.strands.genomic.coordinate)
    cat("\n")
    
    # Count the percentage of this occurence
    case.both.strands.percent <- case.both.strands.genomic.coordinate[, table(sequence) / nrow(env[[genomic.coordinate]]) * 100]
    by.replicate <- sapply(replicate.names, function(rep) {
      case.both.strands.genomic.coordinate[eval(parse(text = rep)) > 0, table(sequence)]
    })
    
    message(paste("Total percentage (out of total unique cases) of such occurence is",
                  sum(case.both.strands.percent), "%"))
    cat("The percentage breakdown is as follow:\n")
    print( case.both.strands.percent )
    
    if(replicate.exist) {
      # Count the percentage of this occurence in each replicate
      cat("\nThe percentage of replicates displaying this behaviour are as follow:\n")
      print(case.both.strands.cnt.by.rep)
      
      ## If there is a count mismatch between complementary sequences, that means the strands are mixed
      ## and matched to produce this observation (same genomic coordinates on both plus and minus strand)
      if ( sum(case.both.strands.cnt.by.rep[c("A"),] != case.both.strands.cnt.by.rep[c("T"),]) > 0 | 
           sum(case.both.strands.cnt.by.rep[c("C"),] != case.both.strands.cnt.by.rep[c("G"),]) > 0 ) {
        
        case.both.strands.cnt.by.rep <- sapply(replicate.names, function(rep) {
          idx <- duplicated(env[[genomic.coordinate]][ !is.na(eval(parse(text = rep))), 1:3]) | 
            duplicated(env[[genomic.coordinate]][ !is.na(eval(parse(text = rep))), 1:3],
                       fromLast = T)
          tab <- env[[genomic.coordinate]][ !is.na(eval(parse(text = rep))) ][idx, table(sequence) /  
                                                                         total.case.reps * 100]
          
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
        cat("Here is the real percentage within each replicate where case occurs at the same position of both strands.\n")
        print(case.both.strands.cnt.by.rep)
        cat("\n")
      }
    }
    
  }
  
  #-------------------------------------------------------------------------------------------------------
  # Venn/Euler Diagram of replicates
  
  if (replicate.exist){
    
    n <- length(replicate.names)
    
    setops <- sapply(1:n, function(i) {
      
      rep.combn <- combn(1:n, i)
      
      rep.combn.names <- apply(rep.combn, 2, paste, collapse="&")
      
      rep.combn.count <- rep(NA, ncol(rep.combn))
      names(rep.combn.count) <- rep.combn.names
      
      for (j in 1:length(rep.combn.names)) {
        
        reps <- rep.combn[, j]
        filter <- paste0("(!is.na(replicate_", reps, "))", collapse = " & ")
        
        rep.combn.count[rep.combn.names[j]] <- env[[genomic.coordinate]][eval(parse(text = filter)), .N]
      }
      
      return(rep.combn.count)
    })
    
    setops <- unlist(setops)
    
    diagram <- venneuler(setops)
    
    png(paste0(output, "/replicates_venneuler.png"))
    plot(diagram)
    dev.off()
    
    # auto plot labeling. hmmm...
    
    # calculate percentage
    setops.percent <- setops / total.union * 100
    
    cat("Damage summary (percentage out of total unique sites)\n")
    print(setops.percent)

  }
  
  
}