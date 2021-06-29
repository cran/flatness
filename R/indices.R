#' Compute flatness indices for (rank) histograms
#'
#' S3 generic function that computes and returns indices of flatness of one or
#' several (rank) histograms, presented in Wilks (2019).
#'
#' Currently the implemented flatness indices are the chi-square statistics, the
#' reliability index and the entropy.
#'
#' @param rkhists an object containing (rank) histograms
#' @param ... other arguments
#' @return the expected returned object is a matrix, with one flatness index in
#' each column for each rank histogram (row-wise). The columns should be named,
#' with "chisq" for the chi-square statistics, "RI" for the reliability index
#' and "entropy" for the entropy.
#' @references
#' Wilks, D. S. "Indices of rank histogram flatness and their sampling
#' properties." \emph{Monthly Weather Review} 147.2 (2019): 763-769.
#' doi:https://doi.org/10.1175/MWR-D-18-0369.1
#' @export
flatness_indices <- function(rkhists, ...) {
  UseMethod("flatness_indices")
}

#' Default method for S3 generic function \code{flatness_indices}
#'
#' Just generate an error.
#'
#' @inheritParams flatness_indices
#' @export
#' @return No return value, called for side effects (generate an error message).
flatness_indices.default <- function(rkhists, ...) {
  stop(paste(c("No methods for function 'flatness_indices' for objects of
               class", class(rkhists))), collapse = " ")
}

#' Method for S3 generic function \code{flatness_indices} and objects with class
#' \code{matrix}
#'
#' This function is the method called when using the S3 generic function
#' \code{flatness_indices} with an object of class \code{matrix}.
#'
#' The input matrices is supposed to contain one rank histograms in each row,
#' unless \code{col_wise == TRUE}.
#'
#' @inheritParams flatness_indices
#' @param col_wise a logical. Are the rank histograms stored column wise, that
#' is one rank histogram in each column?
#' @export
#' @return See \emph{Value} in function \code{flatness_indices}
flatness_indices.matrix <- function(rkhists, col_wise = FALSE, ...) {
  k <- NCOL(rkhists)
  chisq <- apply(rkhists, 1, function(x) {
    n <- sum(x)
    e <- n/k
    sum((x - e)^2/e)
  })
  RI <- apply(rkhists, 1, function(x) {
    n <- sum(x)
    sum(abs(x/n - 1/k))
  })
  entropy <- apply(rkhists, 1, function(x){
    n <- sum(x)
    -1/log(k)*sum(x/n * log(x/n))
  })
  indices <- cbind(chisq = chisq, RI = RI, entropy = entropy)
  return(indices)
}


