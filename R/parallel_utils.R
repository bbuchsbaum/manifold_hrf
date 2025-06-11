# Parallel processing helpers

#' Check if running on Windows
#' @return Logical TRUE if on Windows
#' @keywords internal
.is_windows <- function() {
  .Platform$OS.type == "windows"
}

#' Internal parallel lapply helper
#'
#' Applies a function over a set of elements using base parallel backends when
#' requested. If \code{n_jobs == 1}, a regular \code{lapply} is used. On
#' Unix-like systems, \code{parallel::mclapply} is used for multi-core processing.
#' On Windows, a socket cluster is created for compatibility.
#'
#' @param X Vector or list of elements to iterate over.
#' @param FUN Function to apply to each element.
#' @param n_jobs Number of parallel workers. Set to 1 for sequential execution.
#' @return A list of results from applying \code{FUN}.
#' @keywords internal
.parallel_lapply <- function(X, FUN, n_jobs = 1) {
  if (n_jobs == 1) {
    return(lapply(X, FUN))
  }
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("Package 'parallel' not available. Falling back to sequential execution.")
    return(lapply(X, FUN))
  }
  
  if (.is_windows()) {
    # Use socket cluster on Windows
    cl <- parallel::makeCluster(n_jobs)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, X, FUN)
  } else {
    # Use forking on Unix-like systems
    parallel::mclapply(X, FUN, mc.cores = n_jobs)
  }
}
