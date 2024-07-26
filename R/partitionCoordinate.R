#' Partition overlapping or continuous regions.
#'
#' Table must have start and end columns. The mechanism is similar to the
#'    disjoin function from GenomicRanges but the end coordinate is different.
#'
#' @param coor Coordinate `data.table`.
#' @return Partitioned coordinate `data.table`.
#' 
#' @importFrom data.table setorderv
#' 
#' @export
partitionCoordinate <- function(coor) {

  # Check necessary columns.
  if (!all(c("start", "end") %in% names(coor)))
    stop("Coordinate table must have column start and end.")

  has.col.chr <- "chromosome" %in% names(coor)
  has.col.strand <- "strand" %in% names(coor)
  group.by.cols <- c(if(has.col.chr) "chromosome", if(has.col.strand) "strand")

  # Sort the table
  setorderv(coor, c(group.by.cols, "start", "end"))

  # Locate the overlapping or continuous regions
  coor[, group := cumsum(c(1, cummax(utils::head(end, -1)) - utils::tail(start, -1) < -1)),
       by = eval(group.by.cols)]

  # Make partition of the overlapping regions in each group
  coor2 <- coor[, {
    pos <- union(start, end) |> sort()
    .(start = pos[-length(pos)], end = pos[-1])
  }, by = eval(c(group.by.cols, "group"))]

  # Remove column group
  coor[, group := NULL]
  coor2[, group := NULL]

  return(coor2)
}
