% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rkhist.R
\name{boxplot.rkhist}
\alias{boxplot.rkhist}
\title{Methods for function \code{boxplot} of an S3 object of class \code{rkhist}}
\usage{
\method{boxplot}{rkhist}(x, mini = min(0, min(x)), what = "counts", ...)
}
\arguments{
\item{x}{an \code{rkhist} S3 object containing the rank histograms whose
boxplot must be drawn. See function \code{flatness::rkhist}.}

\item{mini}{a numeric. The minimum value in the \code{ylim} argument
used to plot. Relevant only when plotting counts, ignored otherwise.}

\item{what}{a character string taking its value in
\code{c("counts", "percents", "proportions")}. What must be printed.}

\item{...}{further arguments passed to function
\code{graphics::boxplot.matrix}}
}
\value{
a list as in \code{graphics::boxplot}.
}
\description{
For each rank in an \code{rkhist} object, draws a boxplot of the counts for
all ensembles. A line at the count value for a perfectly flat rank histogram
is also drawn.
}
