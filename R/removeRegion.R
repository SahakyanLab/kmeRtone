#' Remove overlapping region in coordinate `data.table`.
#' 
#' Any "coor" that overlap within the "region" will be removed
#'    e.g. region = 10-20 and coor = 1-30
#'         The results will be: coor = 1-10, 20-30
#'         The coor still overlap one base at the terminal. This is done to
#'         produce exact result as the previous MPhil research.
#' 
#' @param coor Coordinate `data.table`.
#' @param region A `data.table` of region coordinate to be removed.
#' 
#' @return New coordinate `data.table` with the regions removed.
#' 
#' @importFrom data.table rbindlist setkeyv .N
#' @importFrom stringi stri_replace_first_fixed stri_replace_all_regex 
#'  stri_split_fixed stri_sub stri_locate_first_regex stri_replace_first_regex
#'  stri_match_all_regex
#'
#' @export
removeRegion <- function(coor, region) {

  has.col.chr <- "chromosome" %in% names(coor)
  has.col.strand <- "strand" %in% names(coor)
  
  group.by.cols <- c(if(has.col.chr) "chromosome", if(has.col.strand) "strand")
  
  # Combine the coor and the region
  coor <- rbind(coor[, .SD, .SDcols = c(group.by.cols, "start", "end")],
                region[, .SD, .SDcols = c(group.by.cols, "start", "end")])
  setkeyv(coor, c(group.by.cols, "start", "end"))

  # Locate continuous region coordinates and assign group
  coor[, group := cumsum(c(1, cummax(utils::head(end, -1)) - utils::tail(start, -1) < -1)),
       by = eval(group.by.cols)]

  # Make partition of the overlapping regions in each group
  coor <- coor[, {
    pos <- union(start, end) |> sort()
    .(start = pos[-length(pos)], end = pos[-1])
  }, by = eval(c(group.by.cols, "group"))][, group := NULL]

  # Shrink the partitions by one on each of their terminal
  coor[, `:=`(start = start + 1, end = end - 1)]

  # Add original region table to the main coor table again.
  coor <- rbind(coor, region[, .SD, .SDcols = c(group.by.cols, "start", "end")])

  # Find overlaps i.e. shrunken partitions within the case zones
  setkeyv(coor, c(group.by.cols, "start", "end"))
  coor[, group := cumsum(c(1, cummax(utils::head(end, -1)) - utils::tail(start, -1) < 0)),
       by = eval(group.by.cols)]

  # Remove group with more than one partition i.e. the overlaps
  coor <- coor[, if (.N == 1) .(start, end),
               by = eval(c(group.by.cols, "group"))]

  ## Expand the partition back
  coor[, `:=`(start = start - 1, end = end + 1)]

  # Merge continuous partitions
  coor <- mergeCoordinate(coor)

  return(coor)
}
