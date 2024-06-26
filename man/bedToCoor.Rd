% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bedToCoor.R
\name{bedToCoor}
\alias{bedToCoor}
\title{Convert a BED file to chromosome-separated csv files.}
\usage{
bedToCoor(bed.path, output.path = "coordinate/", compress = TRUE)
}
\arguments{
\item{bed.path}{A path to a BED file.}

\item{output.path}{Output directory path. It should be an empty directory.
Default to coordinate/}

\item{compress}{Logical. If TRUE, compress the output CSV files. Default to TRUE.}
}
\value{
None
}
\description{
Convert a BED file to chromosome-separated csv files.
}
