# Logger utilities for M-HRF-LSS

#' Create a simple logger
#'
#' The logger accumulates messages, warnings and decisions during the
#' analysis pipeline. Messages are timestamped for reproducibility.
#'
#' @details
#' **Thread Safety Warning**: This logger is NOT thread-safe. When using 
#' parallel processing (e.g., with future, parallel, or foreach), logger 
#' methods should only be called from the main thread. Calling logger 
#' methods from worker threads may result in lost messages or race conditions.
#' 
#' For parallel processing scenarios, consider:
#' \itemize{
#'   \item Collecting messages in each worker and returning them to the main thread
#'   \item Using the logger only before and after parallel sections
#'   \item Implementing a thread-safe logging solution if concurrent logging is critical
#' }
#'
#' @return An object of class `mhrf_logger` with methods `add` and `get`.
#' @export
create_logger <- function() {
  # Use an environment to properly encapsulate state
  log_env <- new.env(parent = emptyenv())
  log_env$entries <- character()
  
  add <- function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_env$entries <- c(log_env$entries, sprintf("[%s] %s", stamp, msg))
    invisible(NULL)
  }
  
  get <- function() {
    log_env$entries
  }
  
  structure(list(add = add, get = get), class = "mhrf_logger")
}

#' @export
print.mhrf_logger <- function(x, ...) {
  cat("M-HRF-LSS Log\n")
  entries <- x$get()
  if (length(entries) == 0) {
    cat("(empty)\n")
  } else {
    cat(paste(entries, collapse = "\n"), "\n")
  }
  invisible(x)
}
