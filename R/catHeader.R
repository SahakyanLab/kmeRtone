catHeader <- function(msg) {

  space <- rep(" ", (60 - nchar(msg)) / 2)
  cat(rep("-", 60), "\n", sep = "")
  cat(space, msg, space, "\n", sep = "")
  cat(rep("-", 60), "\n", sep = "")

}
