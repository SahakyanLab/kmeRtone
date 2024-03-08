#' Write FASTA files.
#'
#' @param seqs A vector or list of sequences with header name. If it is a list,
#'    it must only contain one single sequence string for every element e.g.
#'    list(chr1 = "NNNNNNNN") not list(chr1 = c("NNNNNN", "AAAAAA"))
#' @param fasta.path A path to a FASTA file.
#' @param append Boolean. Default is FALSE. If TRUE, will append the results to existing file. 
#'
#' @return None
#' 
#' @importFrom data.table fwrite rbindlist
#'
#' @export
writeFASTA <- function(seqs, fasta.path, append = FALSE) {

  seqs.dt <- lapply(seq_along(seqs), function(i) {

    start <- seq(1, nchar(seqs[i]), 60)
    end <- start + 60 - 1
    end[end > nchar(seqs[i])] <- nchar(seqs[i])

    header.name <- names(seqs[i])
    if (!grepl("^>", header.name)) header.name <- paste0(">", header.name)

    seq.dt <- data.table(c(header.name, stri_sub(seqs[i], start, end)))

    return(seq.dt)
  }) |> rbindlist()

  fwrite(seqs.dt, fasta.path, col.names = FALSE, append = append)
}
