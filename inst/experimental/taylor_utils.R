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
  stopifnot(is.numeric(order) && length(order) == 1 && order >= 0)
  stopifnot(is.numeric(x) && length(x) == 1 && is.finite(x))
  stopifnot(is.numeric(h) && length(h) == 1 && h > 0)
  
  if (order == 0) {
    result <- tryCatch({
      f(x)
    }, error = function(e) {
      stop(sprintf("Function evaluation failed at x = %g: %s", x, e$message))
    })
    if (!is.finite(result)) {
      stop(sprintf("Function returned non-finite value at x = %g", x))
    }
    return(result)
  }
  
  g <- function(z) {
    f_plus <- tryCatch({
      f(z + h)
    }, error = function(e) {
      stop(sprintf("Function evaluation failed at x = %g: %s", z + h, e$message))
    })
    
    f_minus <- tryCatch({
      f(z - h)
    }, error = function(e) {
      stop(sprintf("Function evaluation failed at x = %g: %s", z - h, e$message))
    })
    
    if (!is.finite(f_plus) || !is.finite(f_minus)) {
      stop(sprintf("Function returned non-finite values near x = %g", z))
    }
    
    (f_plus - f_minus) / (2 * h)
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
  stopifnot(is.numeric(order) && length(order) == 1 && order >= 0)
  stopifnot(is.numeric(x0) && length(x0) == 1 && is.finite(x0))
  stopifnot(is.numeric(h) && length(h) == 1 && h > 0)
  
  # Check if order is too large (factorial overflow)
  if (order > 170) {
    stop("Order too large, factorial will overflow")
  }
  
  coeffs <- tryCatch({
    sapply(0:order, function(i) {
      deriv <- numeric_derivative(f, x0, i, h)
      fact <- factorial(i)
      if (!is.finite(fact)) {
        stop(sprintf("Factorial overflow at order %d", i))
      }
      coeff <- deriv / fact
      if (!is.finite(coeff)) {
        warning(sprintf("Non-finite coefficient at order %d", i))
      }
      coeff
    })
  }, error = function(e) {
    stop(sprintf("Failed to compute Taylor coefficients: %s", e$message))
  })
  
  return(coeffs)
}

#' Evaluate a Taylor series at given points
#'
#' @param coeff Numeric vector of Taylor coefficients.
#' @param x0 Expansion point used for the coefficients.
#' @param x Numeric vector of evaluation points.
#' @return Numeric vector of the Taylor polynomial evaluated at `x`.
#' @export
evaluate_taylor <- function(coeff, x0, x) {
  # Input validation
  if (!is.numeric(coeff) || length(coeff) == 0) {
    stop("coeff must be a non-empty numeric vector")
  }
  if (!is.numeric(x0) || length(x0) != 1 || !is.finite(x0)) {
    stop("x0 must be a single finite numeric value")
  }
  if (!is.numeric(x) || length(x) == 0) {
    stop("x must be a non-empty numeric vector")
  }
  
  # Check for non-finite coefficients
  if (any(!is.finite(coeff))) {
    warning("Non-finite coefficients detected, replacing with 0")
    coeff[!is.finite(coeff)] <- 0
  }
  
  # Compute powers (0 to length(coeff)-1)
  n <- seq_along(coeff) - 1
  
  # Evaluate polynomial at each point
  result <- tryCatch({
    sapply(x, function(xx) {
      if (!is.finite(xx)) {
        warning(sprintf("Non-finite evaluation point: %s", xx))
        return(NA_real_)
      }
      
      # Compute (x - x0)^n for each power
      x_diff <- xx - x0
      
      # Check for potential overflow
      if (abs(x_diff) > 1 && length(coeff) > 10) {
        # For large x_diff and high order, check for overflow
        max_power <- n[length(n)]
        if (abs(x_diff)^max_power == Inf) {
          warning(sprintf("Polynomial evaluation overflow at x = %g (order %d)", xx, max_power))
          return(NA_real_)
        }
      }
      
      # Compute powers
      powers <- x_diff^n
      
      # Check for non-finite powers
      if (any(!is.finite(powers))) {
        warning(sprintf("Non-finite powers at x = %g", xx))
        powers[!is.finite(powers)] <- 0
      }
      
      # Compute weighted sum
      result <- sum(coeff * powers)
      
      if (!is.finite(result)) {
        warning(sprintf("Non-finite result at x = %g", xx))
        return(NA_real_)
      }
      
      result
    })
  }, error = function(e) {
    stop(sprintf("Failed to evaluate Taylor polynomial: %s", e$message))
  })
  
  return(result)
}
