# Parallel processing helpers

#' Internal parallel lapply helper
#'
#' Applies a function over a set of elements using base parallel backends when
#' requested. If \code{n_jobs == 1}, a regular \code{lapply} is used. On
#' systems where the \pkg{parallel} package is available, \code{parallel::mclapply}
#' is used for multi-core processing.
#'
#' @param X Vector or list of elements to iterate over.
#' @param FUN Function to apply to each element.
#' @param n_jobs Number of parallel workers. Set to 1 for sequential execution.
#' @return A list of results from applying \code{FUN}.
#' @keywords internal
.parallel_lapply <- function(X, FUN, n_jobs = 1) {
  if (n_jobs == 1) {
    lapply(X, FUN)
  } else {
    if (requireNamespace("parallel", quietly = TRUE)) {
      parallel::mclapply(X, FUN, mc.cores = n_jobs)
    } else {
      lapply(X, FUN)
    }
  }
}
