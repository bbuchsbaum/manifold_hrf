# Utility functions

#' Default value operator
#'
#' Returns left hand side if not NULL, otherwise right hand side
#' @param a value to check
#' @param b default value
#' @return a if not NULL else b
#' @export
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

#' Internal parallel lapply helper
#'
#' Applies a function over a set of elements using parallel backends
#' when requested. The function first attempts to use
#' \pkg{future.apply} with a multisession plan, falling back to
#' \code{parallel::mclapply} on Unix-like systems. If neither backend is
#' available or \code{n_jobs = 1}, a regular \code{lapply} is used.
#'
#' @param X Vector or list of elements to iterate over.
#' @param FUN Function to apply to each element.
#' @param n_jobs Number of parallel workers. Set to 1 for sequential
#'   execution.
#' @return A list of results from applying \code{FUN}.
#' @keywords internal
.parallel_lapply <- function(X, FUN, n_jobs = 1) {
  if (n_jobs > 1) {
    if (requireNamespace("future.apply", quietly = TRUE)) {
      oplan <- future::plan(future::multisession, workers = n_jobs)
      on.exit(future::plan(oplan), add = TRUE)
      future.apply::future_lapply(X, FUN)
    } else if (.Platform$OS.type != "windows" &&
               requireNamespace("parallel", quietly = TRUE)) {
      parallel::mclapply(X, FUN, mc.cores = n_jobs)
    } else {
      lapply(X, FUN)
    }
  } else {
    lapply(X, FUN)
  }
}

#' Validate and standardize ridge penalty parameter
#'
#' Ensures a lambda parameter is a non-negative scalar and applies
#' consistent tolerance-based adjustments. Small values below a fixed
#' threshold are coerced to zero with a warning. Unusually large values
#' trigger a warning about potential over-regularization.
#'
#' @param lambda Numeric value provided by the user.
#' @param name Character name of the parameter (for error messages).
#' @return Sanitized lambda value.
#' @keywords internal
.validate_and_standardize_lambda <- function(lambda, name) {
  tol <- 1e-8
  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) || lambda < 0) {
    stop(sprintf("%s must be a non-negative scalar", name))
  }
  if (lambda < tol && lambda > 0) {
    warning(sprintf("%s is near zero (%.2e); treating as 0", name, lambda))
    lambda <- 0
  } else if (lambda > 1) {
    warning(sprintf("%s is large (%.2e); may cause over-regularization", name, lambda))
  }
  lambda
}
