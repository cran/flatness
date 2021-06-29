#' Return a deviate vector with a linear trend
#'
#' Based on the formula given in Jolliffe and Primo (2008) this function returns
#' a vector with the right properties to test the presence of a linear deviate
#' to flatness in a rank histogram.
#'
#' Although the deviate vector is not a rank histogram, this function returns
#' an \code{rkhist} object for the sake of simplicity. This allows for instance
#' to plot the vector.
#'
#' @param k an integer. The number of components.
#' @return An \code{rkhist} object.
#' @references
#' Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' @export
deviate_linear <- function(k) {
  h <- k %/% 2
  res <- k %% 2
  isEven <- as.numeric(res == 0)
  a <- - h - isEven * (h - 1)
  b <- 1 + isEven
  vec <- a + (0:(k-1))*b
  vec <- vec/sqrt(sum(vec^2))
  vec <- matrix(vec, nrow = 1, dimnames = list("linear"))
  class(vec) <- c("rkhist", "matrix")
  return(vec)
}

#' Return a deviate vector with a V-shape trend
#'
#' Based on the formula given in Jolliffe an Primo (2008) this function returns
#' a vector with the right properties to test the presence of a V-shape deviate
#' to flatness in a rank histogram.
#'
#' Although the deviate vector is not a rank histogram, this function returns
#' an \code{rkhist} object for the sake of simplicity. This allows for instance
#' to plot the vector.
#'
#' @references
#' Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' @inheritParams deviate_linear
#' @return An \code{rkhist} object.
#' @export
deviate_V <- function(k) {
  h <- k %/% 2
  res <- k %% 2
  isEven <- res == 0
  if (isEven) {
    a <- h - 1
    b <- 2
    vec1 <- a - (0:(h - 1)) *b
    vec <- c(vec1, rev(vec1))
  } else  {
    a <- h^2
    b <- 2*h + 1
    vec <- a - c(0:(h - 1), (h:0)) * b
  }
  vec <- vec/sqrt(sum(vec^2))
  vec <- matrix(vec, nrow = 1, dimnames = list("V"))
  class(vec) <- c("rkhist", "matrix")
  return(vec)

}

#' Return a deviate vector with a U-shape trend
#'
#' Based on the formula given in Jolliffe and Primo (2008) this function returns
#' a vector with the right properties to test the presence of a U-shape deviate
#' to flatness in a rank histogram.
#'
#' Although the deviate vector is not a rank histogram, this function returns
#' an \code{rkhist} object for the sake of simplicity. This allows for instance
#' to plot the vector.
#'
#' @references
#' Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' @inheritParams deviate_linear
#' @return An \code{rkhist} object.
#' @export
deviate_U <- function(k) {
  h <- k %/% 2
  res <- k %% 2
  isEven <- res == 0
  if (isEven) {
    b <- (k^2 -1)/3
    vec1 <- (k - seq(0, k - 2, by = 2) - 1)^2
    vec <- c(vec1, rev(vec1)) - b
  } else {
    b <- h*(h + 1)/3
    vec1 <- (h - 0:(h - 1))^2
    vec <- c(vec1, round(0, 15), rev(vec1)) - b
  }
  vec <- vec/sqrt(sum(vec^2))
  vec <- matrix(vec, nrow = 1, dimnames = list("U"))
  class(vec) <- c("rkhist", "matrix")
  return(vec)
}

#' Return a deviate vector with deviates at the extreme rank
#'
#' Based on the formula given in Jolliffe and Primo (2008) this function returns
#' a vector with the right properties to test the presence of a deviate to
#' flatness at both extreme ranks in a rank histogram.
#'
#' Although the deviate vector is not a rank histogram, this function returns
#' an \code{rkhist} object for the sake of simplicity. This allows for instance
#' to plot the vector.
#'
#' @references
#' Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' @inheritParams deviate_linear
#' @return An \code{rkhist} object.
#' @export
deviate_ends <- function(k) {
  h <- k %/% 2
  res <- k %% 2
  isEven <- as.numeric(res == 0)
  a <- h - 1 + h * (1 - isEven)
  b <- 1 + (1 - isEven)
  vec <- c(a, rep(-b, k - 2), a)
  vec <- vec/sqrt(sum(vec^2))
  vec <- matrix(vec, nrow = 1, dimnames = list("ends"))
  class(vec) <- c("rkhist", "matrix")
  return(vec)
}

#' Return a deviate vector with a wave-shape trend
#'
#' This function returns a deviate vector with the right properties to test the
#' presence of a deviate to flatness with a sine trend in a rank histogram, as
#' introduced in (Zamo et al. 2021) or (Zamo, 2016).
#'
#' Although the deviate vector is not a rank histogram, this function returns
#' an \code{rkhist} object for the sake of simplicity. This allows for instance
#' to plot the vector.
#'
#' @references
#' \itemize{
#' \item Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' \item Zamo, Michaël, Liliane Bel, and Olivier Mestre. "Sequential aggregation of
#' probabilistic forecasts—application to wind speed ensemble forecasts."
#' Journal of the Royal Statistical Society: Series C (Applied Statistics) 70.1
#' (2021): 202-225. doi:https://doi.org/10.1111/rssc.12455
#' \item Zamo, Michaël. Statistical Post-processing of Deterministic and Ensemble
#' Wind Speed Forecasts on a Grid. Diss. Université Paris-Saclay (ComUE), 2016.}
#' @inheritParams deviate_linear
#' @return An \code{rkhist} object.
#' @export
deviate_wave <- function(k) {
  vec <- sin(2*pi*(0:(k - 1))/(k - 1))
  lin <- deviate_linear(k)
  vec <- vec - sum(vec*lin) *lin
  vec <- vec - sum(vec)
  vec <- vec/sqrt(sum(vec^2))
  vec <- matrix(vec, nrow = 1, dimnames = list("wave"))
  class(vec) <- c("rkhist", "matrix")
  return(vec)
}

#' Check whether a set of vectors obeys the constraints necessary to be used in
#' the Jolliffe-Primo flatness test
#'
#' This function checks whether a set of vectors has the following two
#' properties:
#' \itemize{
#' \item the set is orthonormal
#' \item each vector has components summing to zero (that is, it is a deviation)
#' }
#' @param x an \code{rkhist} or matrix object containing the vectors to check
#' (by row).
#' @param verbose a logical. Should the result of the check be displayed?
#' @param tol a positive numeric. The accepted tolerance for the constraints.
#' @return A list with entries
#' \describe{
#' \item{isOK}{TRUE if the set obeys the constraints, FALSE otherwise}
#' \item{tol}{the tolerance allowed on the constraints}
#' \item{crossprod}{the cross product of the vectors}
#' \item{sums}{the sum of the components of each vector}
#' \item{x}{the checked vectors}}
#' @export
is_JP_ready <- function(x, verbose = TRUE, tol = 0.0001) {
  tol <- abs(tol)
  cprod <- crossprod(t(x))
  rSums <- rowSums(x)
  onorm <- is_orthonormal(x, tol = tol)
  nullSum <- all(rSums >= -tol & rSums <= tol)
  isOK <- onorm & nullSum
  if (verbose) {
    if (! isOK) {
      cat("Vector set is not OK (tolerance is: ", tol, ")\n", sep = "")
      cat("Cross product should be unit diagonal:\n")
      print(cprod)
      cat("Row sums should be null:\n")
      print(rSums)
    } else {
      cat("Vector set is OK (tolerance is: ", tol, ")\n", sep = "")
    }
  }
  return(list(isOK = isOK, tol = tol, crossprod = cprod, sums = rSums, x = x))
}

#' Modify a set of \code{rkhist} objects so that they can be used in the
#' Jolliffe-Primo flatness test
#'
#' The Jolliffe-Primo flatness test requires deviate vectors that are
#' orthonormal and that each has components summing to zero. This function
#' ensures this by using, if it is required, the Grahm-Schmidt method to make
#' the set orthonormal and also makes the sum of each vector's components equal
#' to zero.
#'
#' Note this procedure may greatly change the shape of vectors.
#'
#' @param x an \code{rkhist} object.
#' @param verbose a logical. Should the result of the check be displayed?
#' @return An \code{rkhist} object, containing the modified vectors set obeying
#' the requirements for the Jolliffe-Primo flatness test.
#' @export
make_JP_ready <- function(x, verbose = TRUE) {
  if (! is_JP_ready(x, verbose = verbose)$isOK) {
    names <- row.names(x)
    x <- orthonormalize(x)
    x <- t(apply(x, 1, function(x) x - sum(x)))
    row.names(x) <- names
    class(x) <- c("rkhist", "matrix")
  }
  return(x)
}

#' Return a set of vectors with chosen shapes
#'
#' This function returns an \code{rkhist} object containing vectors with chosen
#' shapes or trends. This is intended to be used to apply the Jolliffe-Primo
#' flatness tests of rank histograms (see Jolliffe and Primo, 2008).
#'
#' The convention is that each row of the \code{rkhist} object contains a
#' vector. It is not required that the set be a basis.
#'
#' For each shape in \code{shapes} this function calls a function named
#' 'deviate_shape' with one argument \code{k}. Some pre-coded functions already
#' exist but the user can easily add its own by following this naming
#' convention. The added function must have only one argument \code{k} and
#' return an \code{rkhist} object. It is advised that the returned deviate
#' vector's components should sum to 0 and have a unit module. But this can be
#' imposed by setting the argument \code{constrain} to \code{TRUE}.
#'
#' If \code{constrain == TRUE} the vector set is modified to have the right
#' properties to be used in the Jolliffe-Primo test, through the Grahm-Schmidt
#' method. It is strongly advised to plot the resulting set with function
#' \code{flatness::plot}, since this transformation may greatly change the
#' shape of the original vectors.
#'
#' @param k an integer. The number of possible ranks.
#' @param shapes a vector of character strings. The required shapes of the
#' vectors.
#' @param constrain a logical. If \code{TRUE} the returned set of vectors
#' is constrained to be orthonormal, with each vector having components
#' summing to 0. This is required to use the vectors in the Jolliffe-Primo
#' flatness test.
#' @return An \code{rkhist} object with each row representing a vector of
#' deviation from flatness.
#' @references
#' \itemize{
#' \item Jolliffe, Ian T., and Cristina Primo. "Evaluating rank histograms using
#' decompositions of the chi-square test statistic." \emph{Monthly Weather
#' Review} 136.6 (2008): 2133-2139. doi:https://doi.org/10.1175/2007MWR2219.1
#' \item Zamo, Michaël, Liliane Bel, and Olivier Mestre. "Sequential aggregation
#' of probabilistic forecasts—application to wind speed ensemble forecasts."
#' Journal of the Royal Statistical Society: Series C (Applied Statistics) 70.1
#' (2021): 202-225. doi:https://doi.org/10.1111/rssc.12455
#' \item Zamo, Michaël. Statistical Post-processing of Deterministic and
#' Ensemble Wind Speed Forecasts on a Grid. Diss. Université Paris-Saclay
#' (ComUE), 2016.}
#' @export
#' @examples
#' deviates <- get_deviates(k = 36, shapes = c("linear", "U", "V", "ends", "wave"))
#' plot(deviates)
#' isJPOK <- is_JP_ready(deviates)
#' JPdeviates <- make_JP_ready(deviates)
#' plot(JPdeviates)
#' JPcheck <- is_JP_ready(JPdeviates)
get_deviates <- function(k, shapes = c("linear", "U", "wave"),
                      constrain = FALSE) {
  deviates <- matrix(0, ncol = k, nrow = length(shapes))
  i <- 0
  for (shape in shapes) {
    args <- list(k = k)
    i <- i + 1
    deviates[i, ] <- do.call(paste0("deviate_", shape), args = list(k = k))
  }
  row.names(deviates) <- shapes
  class(deviates) <- c("rkhist", "matrix")
  if (constrain)
    deviates <- make_JP_ready(deviates)
  return(deviates)
}
