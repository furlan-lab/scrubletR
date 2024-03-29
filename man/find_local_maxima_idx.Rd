% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{find_local_maxima_idx}
\alias{find_local_maxima_idx}
\title{Find local maxima indices in a histogram}
\usage{
find_local_maxima_idx(hist)
}
\arguments{
\item{hist}{Numeric vector representing the histogram}
}
\value{
A list containing indices of local maxima in the histogram
}
\description{
This function finds local maxima indices in a histogram, as scipy's argrelmax
fails on plateaus.
}
\keyword{internal}
