#' Add transparency to color.
#'
#' @param cols Colors in hex format or R color code e.g. "red", "black", etc.
#' @param alpha Alpha value.
#' @return Colors with alpha value in hex format.
#'
#' @export
addAlphaCol <- function(cols, alpha) {
  rgb.vals <- col2rgb(cols) / 255
  transp.cols <- rgb(rgb.vals[1, ], rgb.vals[2, ], rgb.vals[3, ], alpha = alpha)
  return(transp.cols)
}
