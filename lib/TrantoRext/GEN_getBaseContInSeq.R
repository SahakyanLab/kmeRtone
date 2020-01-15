## FUNCTION ####################################################################
# This function takes a parsed sequence (such as after the application of the  #
# fasta parser) and equal sized vectors of starting and ending positions. The  #
# function returns either of the G, C, A, T, N and GC contents of the specified#
# segments, inclusive the starting and ending position bases.                  #
################################################################################
getBaseContInSeq <- function(seq, start.pos, end.pos, content.type="GC"){

  if(length(start.pos)!=length(end.pos)){
    stop("ERROR: getBaseContInSeq: lengths of start.pos != end.pos.")
  }

  if(content.type=="G"){
   cont <- sapply(1:length(start.pos), function(i){
                           sum( seq[ start.pos[i]:end.pos[i] ]=="G" )/
                                         (end.pos[i]-start.pos[i]+1)
                           }, simplify=TRUE, USE.NAMES=FALSE)
  }

  if(content.type=="C"){
   cont <- sapply(1:length(start.pos), function(i){
                           sum( seq[ start.pos[i]:end.pos[i] ]=="C" )/
                                         (end.pos[i]-start.pos[i]+1)
                           }, simplify=TRUE, USE.NAMES=FALSE)
  }

  if(content.type=="A"){
   cont <- sapply(1:length(start.pos), function(i){
                           sum( seq[ start.pos[i]:end.pos[i] ]=="A" )/
                                         (end.pos[i]-start.pos[i]+1)
                           }, simplify=TRUE, USE.NAMES=FALSE)
  }

  if(content.type=="T"){
   cont <- sapply(1:length(start.pos), function(i){
                           sum( seq[ start.pos[i]:end.pos[i] ]=="T" )/
                                         (end.pos[i]-start.pos[i]+1)
                           }, simplify=TRUE, USE.NAMES=FALSE)
  }

  if(content.type=="N"){
   cont <- sapply(1:length(start.pos), function(i){
                           sum( seq[ start.pos[i]:end.pos[i] ]=="N" )/
                                         (end.pos[i]-start.pos[i]+1)
                           }, simplify=TRUE, USE.NAMES=FALSE)
  }

  if(content.type=="GC"){
   cont <- sapply(1:length(start.pos), function(i){
                           sum( seq[ start.pos[i]:end.pos[i] ]=="G" |
                                         seq[ start.pos[i]:end.pos[i] ]=="C" )/
                                         (end.pos[i]-start.pos[i]+1)
                           }, simplify=TRUE, USE.NAMES=FALSE)
  }

  return(cont)

}
## FUNCTION ####################################################################
