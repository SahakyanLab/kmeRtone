#' Mix color
#'
#' This is useful to get overlayed colors.
#'
#' @param cols Colors in hex format or R color code e.g. "red", "black", etc.
#' @param alpha Add alpha transparency value.
#' @return New mixed colors in hex format.
#' 
#' @importFrom grDevices col2rgb rgb
#' 
#' @export
mixColors <- function(cols, alpha) {
  rgb.vals <- col2rgb(cols) / 255
  rgb.vals <- rowSums(rgb.vals) / length(cols)
  mixed.cols <- rgb(rgb.vals[1], rgb.vals[2], rgb.vals[3])
  if (!missing(alpha)) mixed.cols <- addAlphaCol(mixed.cols, alpha)
  return(mixed.cols)
}