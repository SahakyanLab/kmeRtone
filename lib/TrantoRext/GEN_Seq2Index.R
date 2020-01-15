################################################################################
# Recursive function to find the pattern index in the lexicographically sorted 
# strings of the same length as the query.seq, composed of only "A", "C","G" and
# "T". The first entry on the ordered list is assumed to have 0 as its index,
# hence, 1 should be added in case of the desired enumeration starts from 1.
# In case there is a non-ACGT letter, the return will be NA!
################################################################################
Seq2Index <- function(query.seq=NULL){

  if( query.seq[1]=="" ){
    return(0)
  } else {
    nch <- length(query.seq)
    # Last character:
    last.char <- query.seq[nch]
    # Updating the string into the one without the last character:
    if(nch==1){
      query.seq <- ""
    } else {
      query.seq <- query.seq[1:(nch-1)]
    }
    
    return(  4*Seq2Index(query.seq=query.seq) + 
              (match(last.char, c("A","C","G","T"))-1)  )
  
  }
  
}
# Seq2Index(query.seq=c("A","A","C","G","C","G","G","C","T","C"))
################################################################################
