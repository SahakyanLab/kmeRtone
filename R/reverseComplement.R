#' Get reverse complement sequence of DNA
#'
#' @param DNA.sequence DNA sequence can be in a form of character vector or
#'     string. Multiple sequences are accepted.
#' @param form Specify the form: "string" of "vector". Default is "string"
#'
#' @return Reverse complementary sequence
#'
#' @importFrom stringi stri_trans_char stri_reverse
#' 
#' @examples
#' \donttest{
#'    reverseComplement("AAAAA")
#'    reverseComplement(c("AAAAA", "CCCCC"))
#'    reverseComplement(c("A", "A", "A", "A"), form = "vector")
#'  }
#' @export
reverseComplement <- function(DNA.sequence, form="string") {

  if (form == "string") {

    # complement
    # equivalent to base function, chartr but a bit faster
    DNA.sequence <- stri_trans_char(DNA.sequence, "ACGT", "TGCA")

    # reverse
    DNA.sequence <- stri_reverse(DNA.sequence)

  } else if (form == "vector") {

    # locate each nucleotide base in the sequence
    idx.A <- DNA.sequence == "A"
    idx.C <- DNA.sequence == "C"
    idx.G <- DNA.sequence == "G"
    idx.T <- DNA.sequence == "T"

    # complement the sequence
    DNA.sequence[idx.A] <- "T"
    DNA.sequence[idx.C] <- "G"
    DNA.sequence[idx.G] <- "C"
    DNA.sequence[idx.T] <- "A"

    # reverse the sequence
    if (is.null(dim(DNA.sequence))) {
      DNA.sequence <- rev(DNA.sequence)
    } else {
      # if input is a matrix where each column is one DNA
      DNA.sequence <- apply(DNA.sequence, 2, rev)
    }

  }

  return(DNA.sequence)
}
