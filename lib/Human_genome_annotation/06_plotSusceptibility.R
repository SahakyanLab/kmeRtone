plotSusceptibility <- function(susceptibility.count, output){
  
  # A helper function to plot susceptibility in 3 levels: susceptible_kmers, high_risk_kmers, low_risk_kmers
  # These 3 levels are actually column names in susceptibility.count
  plotSusceptibilityLevel <- function(level){
    
    # For title in plot
    if(level == "susceptible_kmers") susceptibility.type <- "dyad content"
    else if(level == "high_risk_kmers") susceptibility.type <- "high-risk susceptibility"
    else if(level == "low_risk_kmers") susceptibility.type <- "low-risk susceptibility"
    
    # Change total_kmers to susceptible_kmers for high and low risk
    if(level %in% c("high_risk_kmers", "low_risk_kmers")) total <- "susceptible_kmers"
    else if(level == "susceptible_kmers") total <- "total_kmers"
    
    element.name <- c("upstream", "UTR5", "CDS", "intron", "UTR3", "downstream", "IGR")
    element.label <- c("upstream", "5'-UTR", "CDS", "intron", "3'-UTR", "downstream", "IGR")
    element.col <- c("blue", "green", "red", "pink", "darkgreen", "darkblue", "grey")
    
    susceptibility.count[, fraction_susceptible_kmers := eval(parse(text = level))/eval(parse(text = total))]
    
    susceptibility.count[eval(parse(text = total)) > 0, {
      
      # Boxplot of susceptibility of transcript elements in sense and antisense strands, separated.
      for(strand_ in unique(strand)){
        
        boxplt <- boxplot(fraction_susceptible_kmers[element=="upstream" & strand==strand_],
                          fraction_susceptible_kmers[element=="UTR5" & strand==strand_], 
                          fraction_susceptible_kmers[element=="CDS" & strand==strand_], 
                          fraction_susceptible_kmers[element=="intron" & strand==strand_],
                          fraction_susceptible_kmers[element=="UTR3" & strand==strand_],
                          fraction_susceptible_kmers[element=="downstream" & strand==strand_],
                          fraction_susceptible_kmers[element=="IGR" & strand==strand_],
                          outline = FALSE, # remove outliers
                          names = element.label,
                          ylab = "Fraction of high-risk susceptible kmers",
                          #ylab = "Dyad content",
                          main = gsub("^.", toupper(substr(susceptibility.type, 1,1)),
                                      paste0(susceptibility.type," of transcript elements and intergenic regions of\n",
                                             strand_, " strand"))
        )
        
        # Label median value
        medians <- boxplt$stats[3, 1:7]
        text(x = seq_along(element.name), y = medians, labels = round(medians, digits = 3), pos = 3, cex = 0.8)
        
        # # Outlier counts
        # outlier.pct <- sapply(seq_along(element.name), function(i) {
        #   top <- sum(fraction_susceptible_kmers[element==element.name[i] & strand==strand_] > boxplt$stats[5, i])
        #   bottom <- sum(fraction_susceptible_kmers[element==element.name[i] & strand==strand_] < boxplt$stats[1, i])
        #   percent <- (top + bottom) / length(fraction_susceptible_kmers[element==element.name[i] & strand==strand_]) * 100
        #   return(percent)
        # })
        # 
        # print(outlier.pct)

        # Plot density
        susceptibility.density <- lapply(element.name, function(e){
          
          d <- density(fraction_susceptible_kmers[element==e & strand==strand_])
          
          # Rescale density
          d$y <- d$y * length(fraction_susceptible_kmers[element==e & strand==strand_])
          
          return(d)
        })
        
        density.lim <- unlist(lapply(susceptibility.density, function(d) c(min(d$y), max(d$y))))
        density.lim <- c(min(density.lim), max(density.lim))

        plot(susceptibility.density[[1]], col=element.col[1],
             xlab = "Fraction of susceptible kmers",
             ylab = "Scaled density", ylim = density.lim, lwd = 2,
             main = paste0("Distribution of susceptible kmers in ", strand_, " strand"))
        for(i in seq_along(element.name)[-1]){
          lines(susceptibility.density[[i]], col = element.col[i], lwd = 2)
        }
        legend("topright", legend = element.label, lty = 1, bty = "n", cex = 0.8,
               col = element.col)
      }
      
      # Boxplot of susceptibility of transcript elements in both strands, combined.
      boxplot(fraction_susceptible_kmers[element=="upstream"],
              fraction_susceptible_kmers[element=="UTR5"], 
              fraction_susceptible_kmers[element=="CDS"], 
              fraction_susceptible_kmers[element=="intron"],
              fraction_susceptible_kmers[element=="UTR3"],
              fraction_susceptible_kmers[element=="downstream"],
              fraction_susceptible_kmers[element=="IGR"],
              outline = FALSE, # remove outliers
              names = element.label,
              #ylab = "Fraction of susceptible kmers",
              ylab = "Dyad content",
              main = gsub("^.", toupper(substr(susceptibility.type, 1,1)),
                          paste0(susceptibility.type," of transcript elements and intergenic regions of\n",
                                 "sense and antisense strand"))
              )
      
      # Label median value
      medians <- sapply(element.name, function(e) median(fraction_susceptible_kmers[element==e]))
      text(x = seq_along(element.name), y = medians, labels = round(medians, digits = 3), pos = 3, cex = 0.8)
      
      strand.comparison <- sapply(element.name, function(element_){
        
        sense <- mean(fraction_susceptible_kmers[element==element_ & strand=="sense"])
        antisense <- mean(fraction_susceptible_kmers[element==element_ & strand=="antisense"])
        
        comparison <- log(sense/antisense, base = 2)
        
        return(comparison)
      })
      
      # Barplot of sense/antisense in log scale
      barplot(strand.comparison,
              ylab = "log2 ratio of sense and antisense high susceptibility",
              #ylab = "log2 ratio of sense and antisense dyad content",
              yaxt = "n", ylim = c(-0.6, 0.1),
              names.arg = element.label,
              main = paste0("Comparison of average ",susceptibility.type," between\nsense and antisense strands"))
      axis(side = 2, at = seq(-0.6, 0.1, 0.1), labels = seq(-0.6, 0.1, 0.1))
      
      NULL
    }]
    
    susceptibility.count[, fraction_susceptible_kmers := NULL]
  }
  
  # Boxplots and barplot
  plotSusceptibilityLevel("susceptible_kmers")
  plotSusceptibilityLevel("high_risk_kmers")
  plotSusceptibilityLevel("low_risk_kmers")
  
  # Correlation plot between all susceptible kmers vs. high/low
  setkey(susceptibility.count, strand, element)
  susceptibility.count[susceptible_kmers > 0, {
    
    plot(susceptible_kmers/total_kmers, high_risk_kmers/susceptible_kmers, main = paste(element, strand),
         #xlab = "Fraction of Susceptible Kmers", ylab = "Fraction of High-risk Kmers")
    xlab = "Dyad content", ylab = "Fraction of High-risk Kmers")
    
  }, by = .(strand, element)]
}