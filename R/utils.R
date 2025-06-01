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

#' Check RAM feasibility for trial precomputation
#'
#' Estimates expected memory usage for storing trial-by-voxel matrices and
#' compares it to a user-provided limit.
#'
#' @param T_trials Number of trials.
#' @param V Number of voxels.
#' @param ram_limit_GB Memory limit in gigabytes.
#'
#' @return Logical indicating whether precomputation is feasible.
#' @keywords internal
check_ram_feasibility <- function(T_trials, V, ram_limit_GB) {
  expected_GB <- (T_trials * V * 8) / 1e9
  feasible <- expected_GB < ram_limit_GB
  if (!feasible) {
    message(sprintf(
      "Estimated memory %.2f GB exceeds limit %.2f GB - using lazy evaluation for trial regressors.",
      expected_GB, ram_limit_GB
    ))
  }
  feasible
}
