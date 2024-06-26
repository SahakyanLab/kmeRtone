% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countChoppedKmers.R
\name{countChoppedKmers}
\alias{countChoppedKmers}
\title{Function chops k-mers within specified ranges of a sequence and counts
them. It uses either a substring method or functionalities from the
Biostrings package.}
\usage{
countChoppedKmers(dna.seq, starts, ends, k, method = "auto")
}
\arguments{
\item{dna.seq}{A string of sequence.}

\item{starts}{Start positions.}

\item{ends}{End positions.}

\item{k}{Size of kmer.}

\item{method}{Method: "Biostrings" or "substring". Default is Biostrings.}
}
\value{
A k-mer-named vector of counts.
}
\description{
Function chops k-mers within specified ranges of a sequence and counts
them. It uses either a substring method or functionalities from the
Biostrings package.
}
