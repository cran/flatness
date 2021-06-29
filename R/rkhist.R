# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Return the rank histogram of an observation in an ensemble
#' forecast
#'
#' This S3 generic function is intended to compute the rank of each observation
#' when pooled with its associated ensemble forecast and return the count in
#' each rank (i.e. the rank histogram).
#'
#' For new methods, the output should be an object of class
#' \code{c("rkhist", "matrix")}, with one rank histogram in each row. Rows may
#' be named.
#'
#' @param fcst an object containing the ensemble forecasts.
#' @param obs an object containing the observation associated to the forecast in
#' \code{fcst}.
#' @param ... additional arguments.
#' @return An S3 object of class \code{rkhist} (indeed a matrix containing the
#' count for each rank, with \code{class} "rkhist"). Each row of the matrix
#' contains the counts for one rank histogram.
#' @examples
#' set.seed(42)
#' N <- 1000
#' M <- 20
#' fcst <- matrix(rnorm(N*M), ncol = M)
#' fcst2 <- matrix(rnorm(N*M), ncol = M)
#' obs <- rnorm(N)
#' ## Computation of one rank histogram
#' # Named
#' rkh <- rkhist(fcst, obs, names = "a")
#' print(rkh)
#' plot(rkh)
#' # Unnamed
#' rkh2 <- rkhist(fcst2, obs)
#' print(rkh2)
#' plot(rkh2)
#'
#' ## Computation of two rank histograms, from a list of forecasts, with the
#' ## same observation vector
#' fcstsl <- list(fcst, fcst2)
#' rkhsl <- rkhist(fcstsl, obs, names = c("a", NA))
#' print(rkhsl)
#' plot(rkhsl)
#'
#' ## Concatenation of two rank histograms, with different names
#' rkhs <- rbind_rkhists(rkh, rkh2, names = letters[3:4])
#' rownames(rkhs)
#' print(rkhs)
#' plot(rkhs)
#' @export
rkhist <- function(fcst, obs, ...) {
  UseMethod("rkhist")
}

#' Default method for S3 generic function \code{rkhist}
#'
#' Just generate an error.
#'
#' @inheritParams rkhist
#' @return No return value, called for side effects (generate an error message).
rkhist.default <- function(fcst, obs, ...) {
  stop(paste(c("No methods for function 'rkhist' for objects of class",
             class(fcst))), collapse = " ")
}

#' Method of S3 generic function \code{rkhist} for \code{matrix} objects
#'
#' This is the method called when the \code{fcst} argument in function
#' \code{rkhist} is a matrix.
#'
#' @param fcst an object containing the ensemble forecasts. It must be a matrix,
#' each row containing a forecast.
#' @param obs an object containing the observation associated to the forecast in
#' \code{fcst}. It can be a vector or a matrix, whose length is the same as
#' the number of row in \code{fcst}
#' @param names a character vector. The row names in the returned \code{rkhist}
#' object.
#' @param ... additional arguments, passed to function \code{base::rank}.
#' @return An S3 object of class \code{rkhist} (indeed just 1-row matrix
#' containing the count for each rank, with \code{class == "rkhist"}).
rkhist.matrix <- function(fcst, obs, names = NULL, ...) {
  obs <- as.vector(obs)
  if (length(obs) != NROW(fcst))
    stop("Unequal number of forecasts and observations in call to function
          rkhist.matrix")
  fcst <- as.matrix(x = fcst)
  K <- ncol(x = fcst) + 1
  obs <- matrix(data = obs, ncol =  1)
  ranks <- apply(X = cbind(obs, fcst), MARGIN = 1, FUN = rank, ...)[1, ]
  rkh <- matrix(data = tabulate(bin = ranks, nbins = K), nrow = 1)
  if (! missing(names))
    row.names(rkh) <- names
  class(rkh) <- c("rkhist", "matrix")
  return(rkh)
}

#' Method of S3 generic function \code{rkhist} for \code{list} objects
#'
#' This is the method called when the \code{fcst} argument in function
#' \code{rkhist} is a list.
#'
#' @param fcst an object containing the ensemble forecasts. It must be a list,
#' each entry containing the forecast of one ensemble. The forecast of each
#' ensemble is a matrix whose rows contain one forecast.
#' @param obs an object containing the observation associated to the forecast in
#' \code{fcst}. It can be a vector or a list, whose length is the same as
#' the number of row in \code{fcst} If a vector, it is used as the observation
#' vector for all the ensemble in \code{fcst}. If it is a list, each entry is
#' associated to the same entry in \code{fcst}.
#' @param names a character vector. The row names in the returned \code{rkhist}
#' object. If missing, the row names are the names of the list \code{fcst}
#' (if any).
#' @param ... additional arguments, passed to function \code{base::rank}.
#' @return An S3 object of class \code{rkhist} (indeed an N-row matrix
#' containing the count for each rank, with \code{class} "rkhist"). N is the
#' length of \code{fcst}, i.e. the number of ensembles.
rkhist.list <- function(fcst, obs, names = NULL, ...) {
  Nfcst <- length(fcst)
  if (! is.list(obs))
    obs <- lapply(1:Nfcst, function(i) obs)
  rkhists <- t(sapply(1:Nfcst, function(i) {
    rkhist(fcst[[i]], obs[[i]], names[i])
  }))
  result <- rbind_rkhists(rkhists)
  if (! missing(names)) {
    row.names(result) <- names
  } else {
    if (! is.null(names(fcst)))
      row.names(result) <- names(fcst)
  }
  return(result)
}

#' Method of S3 generic function \code{rkhist} for \code{data.frame} objects
#'
#' This is the method called when the \code{fcst} argument in function
#' \code{rkhist} is a data frame.
#'
#' @param fcst an object containing the ensemble forecasts. It must be a data
#' frame, each row containing a forecast.
#' @param obs an object containing the observation associated to the forecast in
#' \code{fcst}. It can be a vector or a matrix, whose length is the same as
#' the number of row in \code{fcst}
#' @param names a character vector. The row names in the returned \code{rkhist}
#' object.
#' @param ... additional arguments, passed to function \code{base::rank}.
#' @return An S3 object of class \code{rkhist} (indeed just 1-row matrix
#' containing the count for each rank, with \code{class == "rkhist"}).
rkhist.data.frame <- function(fcst, obs, names = NULL, ...) {
  fcstm <- data.matrix(fcst)
  rkh <- rkhist(fcst = fcstm, obs = obs, ...)
  return(rkh)
}

#' Method for function \code{plot} an S3 object of class \code{rkhist}
#'
#' Plot a rank histogram stored in an object of S3 class \code{rkhist}, with an
#' horizontal dashed line indicating the expected count for a perfectly flat
#' rank histogram.
#'
#' @param x an \code{rkhist} S3 object containing the rank histogram to
#' plot. See function \code{flatness::rkhist}.
#' @param mini a numeric. The minimum value in the \code{ylim} argument
#' used to plot. Relevant only when plotting counts, ignored otherwise.
#' @param what a character string taking its value in
#' \code{c("counts", "percents", "proportions")}. What must be plotted.
#' @param ... other arguments passed to function \code{base::plot.default}.
#' Some arguments are already used in this function, and may not be changed.
#' @return \code{NULL}, invisibly.
plot.rkhist <- function(x, mini = min(0, min(x)), what = "counts", ...) {
#  oldpar <- par(no.readonly = TRUE)
  K <- NCOL(x)
  Nrkh <- NROW(x)
  if (Nrkh != 1) {
#    cols <- floor(sqrt(Nrkh))
#    rows <- ceiling(Nrkh / cols)
#    par(mfrow = c(rows, cols))
  }
  switch(what,
    counts = {
      ylim <- c(mini, max(x)*1.04)
      ylab <- "Count"
    },
    percents = {
      ylim <- c(0, 100)
      ylab <- "Frequency (%)"
      x <- 100 * x/rowSums(x)
    },
    proportions = {
      ylim <- c(0, 1)
      ylab <- "Proportion"
      x <- x/rowSums(x)
    }
  )
  for (n in 1:Nrkh) {
    rkhist1 <- x[n, ]
    if (! is.null(rownames(x)[n]) && ! is.na(rownames(x)[n])) {
      main <- rownames(x)[n]
    } else {
      main <- ""
    }
    plot(x = 1:K, y = rkhist1, xlab = "Rank", ylab = ylab,
                 type = "n", ylim = ylim, xlim = c(0.5, K +0.5),
                 yaxs = "i", xaxs = "i", main = main, ...)
    for (i in 1:K) {
      rect(i - 0.5, 0, i + 0.5, rkhist1[i], col = "darkgrey")
    }
    h <- switch(what,
           counts = sum(rkhist1)/K,
           percents = 100/K,
           proportions = 1/K)
    abline(h = h, lty = 2, lwd = 2, col = "black")
  }
#  par(oldpar)
  invisible(NULL)
}

#' Print an S3 object of class \code{rkhist}
#'
#' Print an S3 object of class \code{rkhist}, using the same layout as
#' \code{base::print.table}.
#'
#' @param x an \code{rkhist} S3 object containing the rank histogram to
#' print. See function \code{flatness::rkhist}.
#' @param what a character string taking its value in
#' \code{c("counts", "percents", "proportions")}. What must be printed.
#' @param ... further arguments, passed to \code{print}.
#' @return No return value, called for side effects.
print.rkhist <- function(x, what = "counts", ...) {
  switch(what,
         percents = {
           x <- 100 * x/rowSums(x)
         },
         proportions = {
           x <- x/rowSums(x)
         },
         counts = {
         }
  )

  if (NROW(x) == 1) {
    trkhist <- as.table(x[1, ])
    dimnames(trkhist) <- list(rank = 1:NCOL(x))
  } else {
    trkhist <- as.table(x)
    dimnames(trkhist) <- list(forecast = rownames(x), rank = 1:NCOL(x))
  }
  print(trkhist, ...)
  invisible(NULL)
}

#' Methods for function \code{boxplot} of an S3 object of class \code{rkhist}
#'
#' For each rank in an \code{rkhist} object, draws a boxplot of the counts for
#' all ensembles. A line at the count value for a perfectly flat rank histogram
#' is also drawn.
#'
#' @param x an \code{rkhist} S3 object containing the rank histograms whose
#' boxplot must be drawn. See function \code{flatness::rkhist}.
#' @param mini a numeric. The minimum value in the \code{ylim} argument
#' used to plot. Relevant only when plotting counts, ignored otherwise.
#' @param what a character string taking its value in
#' \code{c("counts", "percents", "proportions")}. What must be printed.
#' @param ... further arguments passed to function
#' \code{graphics::boxplot.matrix}
#' @return a list as in \code{graphics::boxplot}.
#' @export
boxplot.rkhist <- function(x, mini = min(0, min(x)), what = "counts", ...) {
  K <- NCOL(x)
  switch(what,
         counts = {
           ylim <- c(mini, max(x)*1.04)
           ylab <- "Count"
           h <- sum(x[1, ])/K
         },
         percents = {
           ylim <- c(0, 100)
           ylab <- "Frequency (%)"
           x <- 100 * x/rowSums(x)
           h <- 100/K
         },
         proportions = {
           ylim <- c(0, 1)
           ylab <- "Proportion"
           x <- x/rowSums(x)
           h <- 1/K
         }
  )

  bxp <- boxplot.matrix(x, ylim = ylim, yaxs = "i", xlab = "Rank",
                 ylab = ylab, ...)
  abline(h = h, lwd = 2, lty = 2, col = "black")
  return(bxp)
}

#' Stack objects of class \code{rkhist}
#'
#' Use \code{rbind} to stack several rank histograms (stored in S3 objects of
#' class \code{rkhist}, or matrix, or vectors).
#'
#' @param ... matrix, rkhist objects or vectors to bind by row.
#' @param names a character vector of names to give to each rank histogram
#' (i.e. each row of the resulting matrix). Although its length must be equal to
#' the number of row in the resulting matrix, \code{NA}s are allowed.
#' @return An object of S3 class \code{rkhist}
#' @export
rbind_rkhists <- function(..., names = NULL) {
  rkhists <- rbind(...)
  if (! missing(names))
    rownames(rkhists) = names
  class(rkhists) <- c("rkhist", "matrix")
  return(rkhists)
}
