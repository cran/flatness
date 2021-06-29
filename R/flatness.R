#' flatness: a package to assess the flatness of (rank) histograms
#'
#' The \code{flatness} package offers tools (scores, tests, ...) to  compute
#' histograms and assess whether they are flat.
#'
#' The S3 generic function \code{rkhist} allows to compute one or several rank
#' histograms from ensemble forecasts and corresponding scalar observations. (In
#' Meteorology an ensemble forecast is a set of forecasts for the same variable,
#' aimed at assessing the forecasting uncertainty). It creates an object with
#' class \code{rkhist} that can then be \code{plot}ted and \code{print}ed.
#'
#' Flatness of (rank) histograms may then be tested with function \code{JP_test}
#' that implements the Jolliff-Primo flatness tests. This test requires a
#' set of deviance vectors, some of which can be provided with functions named
#' \code{deviate_XXX}. The user can easily implement its own deviate-returning
#' function (please see details in \code{get_deviates} on how to do this).
#' Functions \code{is_JP_ready} and \code{make_JP_ready} are provided to ensure
#' that a set of deviate vectors meet the requirements to be used in the
#' Jolliffe-Primo tests. The result of the test is stored in an object with
#' class \code{JPtest} that can be \code{print}ed or drawn with the function
#' \code{lattice::levelplot}.
#'
#' Flatness indices can be computed with the S3 generic function
#' \code{flatness_indices}.
#'
#' See the vignette for further details and an illustration with the datasets
#' \code{ensembles} and \code{ppensembles} provided with this package.
#'
#' @docType package
#' @name flatness
NULL
