% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{orthonormalize}
\alias{orthonormalize}
\title{Make a vector set orthonormal}
\usage{
orthonormalize(set)
}
\arguments{
\item{set}{a matrix. The convention used here is that each row of \code{set}
contains a vector.}
}
\value{
A matrix with the same dimension as \code{set} containing an
orthonormal set of vector. The vector are stored in each row.
}
\description{
This function uses the Grahm-Schmidt method to make a set of vector
orthonormal.
}
