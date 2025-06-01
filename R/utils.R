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

#' Adjust HRF Vector for Data Bounds
#'
#' Truncates or pads an HRF vector so that its length does not exceed the
#' available number of time points.
#'
#' @param hrf Numeric vector representing the HRF shape.
#' @param max_timepoints Maximum allowed length.
#' @return HRF vector of length \code{max_timepoints}.
#' @export
adjust_hrf_for_bounds <- function(hrf, max_timepoints) {
  if (!is.numeric(hrf)) {
    stop("hrf must be numeric")
  }

  if (length(hrf) > max_timepoints) {
    warning("HRF truncated to fit within available timepoints")
    hrf[seq_len(max_timepoints)]
  } else if (length(hrf) < max_timepoints) {
    c(hrf, rep(0, max_timepoints - length(hrf)))
  } else {
    hrf
  }
}
