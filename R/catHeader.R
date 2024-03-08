#' Function prints a given message in a formatted header with borders.
#'
#' @param msg message to be printed within the header.
#'
#' @export
catHeader <- function(msg) {

  space <- rep(" ", (60 - nchar(msg)) / 2)
  cat(rep("-", 60), "\n", sep = "")
  cat(space, msg, space, "\n", sep = "")
  cat(rep("-", 60), "\n", sep = "")

}
