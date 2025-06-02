# Logger utilities for M-HRF-LSS

#' Create a simple logger
#'
#' The logger accumulates messages, warnings and decisions during the
#' analysis pipeline. Messages are timestamped for reproducibility.
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
