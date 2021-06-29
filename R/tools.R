#' Make a vector set orthonormal
#'
#' This function uses the Grahm-Schmidt method to make a set of vector
#' orthonormal.
#'
#' @param set a matrix. The convention used here is that each row of \code{set}
#' contains a vector.
#' @return A matrix with the same dimension as \code{set} containing an
#' orthonormal set of vector. The vector are stored in each row.
#' @export
orthonormalize <- function(set) {
  if (is_orthonormal(set))
    return(set)
  if (NROW(set) > 1L) {
    otset <- t(sapply(2:NROW(set), function(i) {
      residuals(lm(t(set[i, , drop = FALSE]) ~ t(set[1:(i-1), , drop = FALSE])))
    }))
    otset <- rbind(set[1, ], otset)
  } else {
    otset <- set
  }
  otset <- t(apply(otset, 1, function(x) x/sqrt(sum(x^2))))
  return(otset)
}

#' Check whether a vector set is orthonormal
#'
#' @inheritParams orthonormalize
#' @param tol a positive numeric. The accepted tolerance for the conditions of
#' orthonormality.
#' @return TRUE if the set of vectors is orthonormal, FALSE otherwise.
#' @export
is_orthonormal <- function(set, tol = 0.0001) {
  tol <- abs(tol)
  crsprd <- crossprod(t(set))
  mdiag <- diag(nrow = NROW(crsprd), ncol = NROW(crsprd))
  diff <- abs(crsprd - mdiag)
  return(all(diff <= tol))
}

#' Multiple statistical hypothesis testing with the Benjamini-Hochberg procedure
#'
#' This function applies the procedure described in Benjamini & Hochberg (1995)
#' for controlling the False Discovery Rate in multiple statistical hypothesis
#' testing.
#'
#' The procedure works as follows. Let N p-values \eqn{p_i}
#' (with \eqn{i=1,\ldots,N}) and a significance level \eqn{\alpha}. The decision
#' threshold is \eqn{p_k} where \eqn{k = max(i;p_{(i)} <= \alpha \frac{i}{N})}
#' where \eqn{p_{(i)}} are the sorted p-values \eqn{p_i}. For every test with
#' \eqn{p_i \le p_k}, the null hypothesis is rejected. By convention
#' \eqn{p_{(0)} = 0}.
#'
#' @param pvalues a vector of p-values, with length N
#' @param alpha a numeric between 0 and 1
#' @param ... same arguments as in function base::sort, except "index.return"
#' which is set to TRUE.
#'
#' @return A named listed with entries:
#' \itemize{
#' \item H0: a logical vector of length N. The nth entry has value TRUE if the
#' null hypothesis associated with the nth p-value is not rejected and FALSE
#' otherwise.
#' \item pk: the decision threshold. A p-value under this threshold leads to
#' rejection of the associated null hypothesis.
#' \item alpha: the chosen significance level
#' \item pvalues: the vector of p-values provided to apply the
#' Benjamini-Hochberg procedure.
#' }
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). "Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing." \emph{Journal of the
#' Royal statistical society: series B (Methodological)}, 57(1), 289-300.
#' doi:https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
#' @export
BH_procedure <- function(pvalues, alpha = 0.01, ...) {
  if ( alpha < 0 || alpha > 1 )
    stop("alpha must be in [0;1]")
  N <- length(pvalues)
  results <- list(H0 = rep(TRUE, N),
                  pk = 0,
                  alpha = alpha,
                  pvalues = pvalues
  )
  pord <- sort(pvalues, ..., index.return = TRUE)
  i <- 1:N
  k <- max(which(pord$x <= i*alpha/N))
  if (! is.infinite(k)) {
    pk <- results$pk <- pord$x[k]
    results$H0[which(pvalues <= pk)] <- FALSE
  }

  return(results)
}
