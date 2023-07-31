#' A R6 class wrapper for data.table
#' 
#' A way to grow data.table in different environment but still retaining access
#'      to it. A temporary fix until data.table developer develop update row
#'      by reference.
#' 
Kmer_Table <- R6::R6Class(
  
  classname = "Kmer_Table",
  public = list(
    
    #' @field DT data.table of k-mers
    DT = data.table(),
    
    #' @description
    #' initialize empty data.table of k-mers
    initialize = function() {
      
      self$DT <- data.table(kmer = character(),
                            pos_strand = numeric(),
                            neg_strand = numeric(),
                            key = "kmer")
      
    },
    
    #' @description
    #' Print method.
    print = function() {
      
      cat("This is R6 wrapper for data.table\n\n")
      print(self$DT)
      
    },
    
    #' @description
    #' Remove unknown base N.
    remove_N = function() {
      
      self$DT <- self$DT[stri_detect_fixed(kmer, "N", negate = TRUE)]
      setkey(self$DT, kmer)
      
    },
    
    #' @description
    #' Filter out k-mers without defined central patterns.
    #' @param central.pattern Central pattern.
    #' @param k Length of k-mer.
    #' @return None.
    filter_central_pattern = function(central.pattern, k) {
      
      self$remove_N()
      
      flank <- (k - unique(stri_length(central.pattern))) / 2
      
      self$DT <- self$DT[
        stri_sub(kmer, flank + 1, length = stri_length(central.pattern)) %in%
          central.pattern
      ]
      
      setkey(self$DT, kmer)
      
    },
    
    #' @description
    #' Update count for existed k-mers in the table.
    #' @param kmers K-mer table with new count to be added to the main table.
    #' @param is.strand.sensitive Does strand polarity matter?
    #' @param strand If yes, what is the strand refers to? "+" or "-".
    #' @return None.
    update_count = function(kmers, is.strand.sensitive, strand) {
      
      # Update kmer.table for existing k-mers in the table
      if (is.strand.sensitive && strand == "-") {
        self$DT[kmers, neg_strand := neg_strand + N]
      } else {
        self$DT[kmers, pos_strand := pos_strand + N]
      }
      
    },
    
    #' @description
    #' Add new rows for new k-mers with their respective counts that is not
    #'      existed yet in the main table.
    #' @param kmers K-mer table with new k-mers to be added to the main table.
    #' @param is.strand.sensitive Does strand polarity matter?
    #' @param strand If yes, what is the strand refers to? "+" or "-".
    #' @return None.
    update_row = function(kmers, is.strand.sensitive, strand) {
      
      # Add additional column to kmers
      if (is.strand.sensitive && strand == "-") {
        setnames(kmers, "N", "neg_strand")
        kmers[, pos_strand := 0]
      } else {
        setnames(kmers, "N", "pos_strand")
        kmers[, neg_strand := 0]
      }
      
      # Bind new k-mers found
      self$DT <- rbind(self$DT, kmers[!self$DT])[
        ,
        `:=`(pos_strand = sum(pos_strand),
             neg_strand = sum(neg_strand)),
        
        by = kmer
      ]
      setkey(self$DT, kmer)
      
    }
    
  )
)