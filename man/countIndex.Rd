% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{countIndex}
\alias{countIndex}
\title{countIndex}
\usage{
countIndex(yIndex, numberBin)
}
\arguments{
\item{yIndex}{vector of indices representing which bins the responses fall into.}

\item{numberBin}{the number of bins for response discretization.}
}
\value{
The function returns a vector indicating the number of responses in each bin, of length \code{numberBin}.
}
\description{
This function counts samples in each bin.
}
