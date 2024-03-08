#' Merge overlapping or continuous regions.
#'
#' Table must have start and end columns. The output is exactly similar to the
#'    reduce function from GenomicRanges.
#'
#' @param coor Coordinate `data.table`.
#' 
#' @return Merged coordinate `data.table`.
#' 
#' @importFrom data.table setorderv
#'
#' @export
mergeCoordinate <- function(coor) {

  # Check necessary columns.
  if (!all(c("start", "end") %in% names(coor)))
    stop("Coordinate table must have column start and end.")

  has.col.chr <- "chromosome" %in% names(coor)
  has.col.strand <- "strand" %in% names(coor)
  group.by.cols <- c(if(has.col.chr) "chromosome", if(has.col.strand) "strand")

  # Sort the table
  setorderv(coor, c(group.by.cols, "start", "end"))

  # Locate the overlapping or continuous regions
  coor[, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < -1)),
       by = eval(group.by.cols)]

  # Merge
  coor2 <- coor[, .(start = min(start), end = max(end)),
                by = eval(c(group.by.cols, "group"))]

  # Remove column group
  coor[, group := NULL]
  coor2[, group := NULL]

  return(coor2)
}
