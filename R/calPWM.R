#' Calculate position weight matrix of overlapping sequences.
#'
#' Simulation of human population is based on single nucleotide variantion.
#'
#' @param kmers A vector of k-mers to overlap.
#' @param pseudo.num Pseudo-number to avoid numerical instability due to lack of
#'     base at a position. Default is zero i.e. no pseudo-number.
#' @param bg.prop Background proportion of bases. Default is c(a = 0.295,
#'     c = 0.205, g = 0.205, t = 0.295) which is observed in human genome.
#' @param output Output matrix type. Options are PCM, PPM, and PWM which refer
#'     to position count/probability/weight matrix. Default is PWM.
#'
#' @return A position count/probability/weight matrix.
#'
#' @export
calPWM <- function(kmers, pseudo.num=0,
  bg.prop=c(a = 0.295, c = 0.205, g = 0.205, t = 0.295), output="PWM") {

  if (length(bg.prop) != 4) {
    stop("Please give proportion for the four DNA bases.")
  } else if (sum(bg.prop) != 1) {
    stop("Total proportion is not equal to one.")
  }

  dt <- stri_split_boundaries(kmers, type = "character", simplify = TRUE) |>
    as.data.table()

  for (j in 1:ncol(dt))
    set(dt, j = j, value = factor(dt[[j]], levels = c("A", "C", "G", "T")))

  # Position count matrix
  pcm <- dt[, sapply(.SD, table)]

  if (stri_trans_tolower(output) == "pcm") return(pcm)

  # Position probability matrix with pseudo number
  ppm <- (pcm + (pseudo.num / 4)) / (length(kmers) + pseudo.num)

  if (stri_trans_tolower(output) == "ppm") {
    attr(ppm, "pseudo.num") <- pseudo.num
    return(ppm)
  }

  # Position weight matrix
  pwm <- log2(ppm / bg.prop[c("a", "c", "g", "t")])

  attr(pwm, "pseudo.num") <- pseudo.num
  attr(pwm, "bg.prop") <- bg.prop

  return(pwm)
}
