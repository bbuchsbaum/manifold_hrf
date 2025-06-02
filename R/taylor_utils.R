# Taylor Series Utility Functions
#
#' Numeric derivative using central differences
#'
#' @param f Function to differentiate.
#' @param x Point at which to compute the derivative.
#' @param order Derivative order (non-negative integer).
#' @param h Step size for finite differences.
#' @return Numeric derivative approximation.
#' @keywords internal
numeric_derivative <- function(f, x, order = 1, h = 1e-4) {
  stopifnot(is.function(f))
  stopifnot(order >= 0)
  if (order == 0) {
    return(f(x))
  }
  g <- function(z) {
    (f(z + h) - f(z - h)) / (2 * h)
  }
  if (order == 1) {
    g(x)
  } else {
    numeric_derivative(g, x, order - 1, h)
  }
}

#' Compute Taylor series coefficients numerically
#'
#' @param f Function to expand.
#' @param x0 Point around which to expand.
#' @param order Maximum order of the series.
#' @param h Step size for finite differences.
#' @return Numeric vector of coefficients of length `order + 1`.
#' @export
numeric_taylor_coefficients <- function(f, x0, order, h = 1e-4) {
  stopifnot(is.function(f))
  stopifnot(order >= 0)
  sapply(0:order, function(i) numeric_derivative(f, x0, i, h) / factorial(i))
}

#' Evaluate a Taylor series at given points
#'
#' @param coeff Numeric vector of Taylor coefficients.
#' @param x0 Expansion point used for the coefficients.
#' @param x Numeric vector of evaluation points.
#' @return Numeric vector of the Taylor polynomial evaluated at `x`.
#' @export
evaluate_taylor <- function(coeff, x0, x) {
  n <- seq_along(coeff) - 1
  sapply(x, function(xx) sum(coeff * (xx - x0) ^ n))
}
