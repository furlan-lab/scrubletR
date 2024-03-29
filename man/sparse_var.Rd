% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{sparse_var}
\alias{sparse_var}
\title{Calculate variance across the specified axis for a sparse matrix}
\usage{
sparse_var(E, axis = 1)
}
\arguments{
\item{E}{A sparse matrix where rows represent observations and columns represent features.}

\item{axis}{An integer specifying the axis along which the variance is calculated (1 for rows, 2 for columns).}
}
\value{
A numeric vector containing the variance values.
}
\description{
This function computes the variance across the specified axis for a sparse matrix.
}
\keyword{internal}
