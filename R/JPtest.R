#' Implementation of the Jolliffe-Primo flatness tests for rank histograms
#'
#' Given a matrix of rank histograms and an orthonormal set of deviate
#' vector(s), this function computes the projection components, test statistics
#' and p-values of the Jolliffe-Primo flatness test for each inputted rank
#' histogram. See Jolliffe and Primo (2008) for details of the method.
#'
#' Note that the test statistics and p-values of the projections over
#' the residual vector (after removing all the projection on the deviates) are
#' also computed and returned.
#'
#' @param rkhists an rkhist, a matrix, or any other object that can be coerced to
#' a matrix. It contains the rank histogram(s) whose flatness must be tested
#' (one in each row).
#' @param deviates the matrix containing the deviate vectors used for testing.
#' Each row contains a deviate vector: the vector set must be orthonormal, and
#' each deviate vector must have its components summing to zero.
#' @param ... further arguments (currently not used).
#' @return  A list (with additional first class \code{JPtest}) with the
#' following entries:
#' \describe{
#' \item{test}{an array containing the result of the Jolliffe-Primo test(s). The
#' first dimension is of length three (the projection over the deviate vectors,
#' the test statistics and the p-values), the second and third dimensions
#' correspond to the rank histogram and the test, respectively}
#' \item{deviates}{the set of deviate vectors used in the test}
#' \item{rkhist}{the tested rank histogram(s) (an \code{rkhist} object).}
#' }
#' @examples
#' require(lattice)
#' require(xtable)
#' M <- 15
#' N <- 100
#' n <- 20
#' fcsts <- vector("list", n)
#' names(fcsts) <- letters[1:n]
#' obs <- rnorm(N)
#' for (i in 1:n) {
#'   fcsts[[i]] <- matrix(rnorm(M*N), ncol = M)
#' }
#' rkhsts <- rkhist(fcsts, obs)
#' deviates <- get_deviates(M + 1)
#' test <- JP_test(rkhsts, deviates)
#' print(test)
#' for (what in c("projections", "statistics", "pvalues")){
#'   levelplot(test, what = what, main = what, rotate = what == "pvalues")
#' }
#' xtable(test$test["pvalues", ,])
#' xtable(t(test$test["pvalues", ,]))
#' @references
#' Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' @export
JP_test <- function(rkhists, deviates, ...) {
  k <- NCOL(deviates) # number of possible ranks
  N <- NROW(deviates) # number of tested shapes
  E <- NROW(rkhists) # number of forecasts

  projections <- tcrossprod(rkhists, deviates)

  residualv <- t(sapply(1:E, function(e) {
    projm <- projections[e, ]
    proj <- colSums(projm * deviates)
    residual <- rkhists[e, ] - proj
    return(residual)
  }))

  statistics <- cbind(projections^2, rowSums(residualv^2))
  residpval <- 1 - pchisq(rowSums(residualv^2), df = k - N - 1)
  pvalues <- 1 - apply(statistics[, - (N + 1), drop = FALSE], 1:2, pchisq,
                       df = 1)
  pvalues <- cbind(pvalues, residpval)

  test <- array(data = NA, dim = c(3, E, N + 1))
  dnames <- list(result = c("projections", "statistics", "pvalues"),
                 rkhist = row.names(rkhists),
                 test = c(vector("character", N), "resid."))
  if (! all(is.null(row.names(deviates))))
    dnames$test[1:N] <- row.names(deviates)
  dimnames(test) <- dnames

  test[1, , 1:N] <- projections
  test[2, , ] <- statistics
  test[3, , ] <- pvalues
  result <- list(deviates = list(deviates),
                 rkhists = list(rkhists),
                 test = test)
  class(result) <- c("JPtest", class(result))
  return(result)
}

#' Print method for a \code{JPtest} object
#'
#' @param x the \code{JPtest} object that must be printed.
#' @param what a character vector. What must be printed? Can contain any or
#' several among c("projections", "statistics", "pvalues").
#' @param ... further arguments. Passed to \code{print} method for
#' matrix objects.
#' @export
#' @return No return value, called for side effects.
print.JPtest <- function(x, what = c("projections", "statistics", "pvalues"),
                       ...) {
  cat("Results of the Jolliffe-Primo test:\n")
  for (w in what) {
    cat(paste0(w, ":\n"))
    print(x$test[w, , ], ...)
  }
  invisible(NULL)
}

#' Levelplot method for data in a \code{JPtest} object
#'
#' Plot a chosen result matrix contained in the \emph{test} entry of a
#' \code{JPtest} object. The underlying function is the
#' \code{lattice::levelplot.matrix} function.
#'
#' @param JPobj the \code{JPtest} object to plot.
#' @param what a character or an integer. Which component of the \code{test}
#' entry in \code{JPobj} should be plotted? Can be any of \emph{projections},
#' \emph{statistics}, \emph{pvalues} or, respectively, 1, 2 or 3.
#' @param rotate a logical. Should the matrix containing the data be transposed
#' before plotting? By default in \code{lattice::levelplot.matrix}, the rows of
#' the plotted matrix correspond to the x-axis.
#' @param plot a logical. Should the lattice object be plotted?
#' @param ... further arguments passed to the \code{lattice::levelplot}
#' function.
#' @return An object of class \code{trellis}, invisibly.
#' @export
levelplot.JPtest <- function(JPobj, what = "pvalues", rotate = TRUE,
                             plot = TRUE, ...)
{
  args <- list(...)
  what <- what[1]
  x <- as.matrix(JPobj$test[what, , ])
  if (rotate) {
    n <- NROW(x)
    x <- t(x)[, n:1, drop = FALSE]
  }
  args$x <- x
  ploplo <- do.call(what = lattice::levelplot, args = args)
  if (plot) {
    print(ploplo, )
  }
  invisible(ploplo)
}

#' #' Method of S3 generic function \code{c} for \code{JPtest} objects
#' #'
#' #' @param ...
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' c.JPtest <- function(...) {
#'   args <- list(...)
#'   JPtests <- args[[which(! names(args) %in% formals(abind::abind))]]
#'   abind_args <- args[[which(names(args) %in% formals(abind::abind))]]
#'   tests <- lapply(JPtests, function(x) x$test)
#'   result <- do.call(abind::abind, args = c(tests, along = 2, abind_args))
#'   return(result)
#' }

#' JPtest <- function(rankhists, basis, ...) {
#'   UseMethod("JPtest")
#' }
#'
#' JPtest.default <- function(rankhists, basis, ...) {
#'   stop(paste("No function JPtest for class", class(rankhists)))
#' }
#'
#' JPtest.rkhist <- function(rankhists, basis, ...) {
#'   rkhv <- as.vector(rankhists)
#'   results <- JPtest.vector(rankhists = rkhv, basis = basis, ...)
#'   return(results)
#' }
#'
#' Jptest.matrix <- function(rankhists, basis, orthonormalize = FALSE) {
#'   results <- apply(rankhists, 1, JPtest, basis = basis, ...)
#'   return(results)
#' }
#'
#' JPtest.vector <- function(rankhists, basis, orthonormalize = FALSE, ...) {
#'   K <- length(rankhists)
#'   n <- sum(rankhists)
#'   e <- rep(n/K, K)
#'   sqe <- sqrt(e)
#'   x <- matrix((rankhists - e)/sqe, ncol = 1)
#'
#'   if (orthonormalize) {
#'     basis <- orthonormalize(basis)
#'   }
#'
#'   u <- basis %*% x
#'   u2 <- u^2
#'   pval <- pchisq(u2, df = 1, lower.tail = FALSE)
#'
#'   results <- cbind(statistics = u2, pvalue = pval)
#'   return(results)
#' }
