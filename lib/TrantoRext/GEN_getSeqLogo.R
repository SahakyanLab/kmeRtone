################################################################################
# This function uses the seqLogo package of R, and generates a seqLogo from a  #
# provided vector of sequences. The function also returns the constructed pwm  #
# matrix. The arguments of the function are as follows:                        #
#                                                                              #
# seq - vector of strings                                                      #
# A vector of strings with the constituent values being the sequences of the   #
# same length.                                                                 #
#                                                                              #
# ic.sclae - a logical value (default = FALSE)                                 #
# An argument, to be passed to the seqLogo function, that defines whether the  #
# logos should be scaled by the information content.                           #
#                                                                              #
# xfontsize and yfontsize (numeric values, default = 11)                       #
# The font sizes to be used in the seqLogo plots.                              #
#                                                                              #
# indices - a numeric vector (default = NULL)                                  #
# If not null, the subset or rearrangement of the bases in the provided se-    #
# quence set will be used for the seqLogo creation.                            #
################################################################################
getSeqLogo <- function(seq, ic.scale=FALSE, xfontsize=11, yfontsize=11, indices=NULL){

  # source("https://bioconductor.org/biocLite.R")
  # biocLite("seqLogo")
  library("seqLogo")

  # determining the lengths of the sequences.
  lengths <- as.vector(unique(nchar(seq)))
  num.of.seq <- length(seq)
  
  if(length(lengths)!=1){ stop("getSeqLogo: the sequence lengths are not equal.") }
  
  if( !is.null(indices[1]) ){ # Indices are not NULL, hence only a subset or a 
                              # rearrangement of the characters in each sequence 
                              # will be examined.
    if( length(indices)>lengths ){stop("getSeqLogo: the index length > sequence lengths.")}
  } else {
    indices <- 1:lengths
  }
  
  prob.matrix <- matrix( 0, nrow=4, ncol=length(indices) )
  dimnames(prob.matrix)[[1]] <- c("A","C","G","T")
  
  for(i in 1:length(indices)){
    chars <- substr(start=indices[i], stop=indices[i], x=seq)
    for(j in c("A","C","G","T")){
      prob.matrix[j,i] <- length(which(chars==j))/num.of.seq
    }
  }

  pwm <- makePWM(prob.matrix)
  #slotNames(pwm)
  seqLogo(pwm, ic.scale=ic.scale, xfontsize=xfontsize, yfontsize=yfontsize)
  return(pwm)
  
}
################################################################################

