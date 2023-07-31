#' Trim out-of-bound coordinates
#' 
#' It operates in two mode: coordinate table with and without chromosome. The
#'    former require `Genome` for the chromosomal sequence length.
#'
#' @field coor Coordinate `data.table`.
#' @field seq.len Sequence length to trim end position.
#' @field genome `Genome` class object.
#' @return Trimmed coordinate `data.table`.
#'
#' @export
trimCoordinate <- function(coor, seq.len, genome) {

  has.col.end <- "end" %in% names(coor)
  has.col.chr <- "chromosome" %in% names(coor)

  # Remove out-of-range coordinates.
  if (has.col.end && coor[end < 1, .N] > 0) coor <- coor[end > 0]
  
  # Remove out-of-range coordinates.
  idx <- coor[, start > if(!has.col.chr) seq.len
                        else genome$get_length(chromosome)]
  if (any(idx)) coor <- coor[!idx]
  
  # Trim overflow start positions.
  coor[start < 1, start := 1]

  # Trim overflow end positions.
  if (has.col.end)
    coor[end > if(!has.col.chr) seq.len else genome$get_length(chromosome),
         end := if(!has.col.chr) seq.len else genome$get_length(chromosome)]

  return(coor)
}