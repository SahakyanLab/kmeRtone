% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initKmerTable.R
\name{initKmerTable}
\alias{initKmerTable}
\title{Initialise k-mer table with all possible k-mers}
\usage{
initKmerTable(k, central.pattern = NULL, split.kmer = FALSE)
}
\arguments{
\item{k}{K-mer size. Limit to 15 because vector size is limited to
.Machine$integer.max. For 9- to 15-mer, the kmer sequence is separated to
two columns (kmer_part1 and kmer_part2) to reduce memory significantly.}

\item{central.pattern}{Central pattern(s) of the k-mer. Default is NULL.}

\item{split.kmer}{Whether to split the k-mer sequence into two parts for large k values. Default is FALSE.}
}
\value{
data.table with 3 columns: kmer, pos_strand, neg_strand
}
\description{
Initialise k-mer table with the following columns: kmer, pos_strand, and
neg_strand
}
