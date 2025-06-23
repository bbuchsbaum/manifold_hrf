#' Null coalescing operator
#'
#' Returns the left-hand side if it is not NULL, otherwise returns the right-hand side.
#'
#' @name grapes-or-or-grapes
#' @param x Left-hand side value
#' @param y Right-hand side value (default)
#' @return \code{x} if not NULL, otherwise \code{y}
#' @examples
#' # Returns 5
#' NULL %||% 5
#' 
#' # Returns 10
#' 10 %||% 5
#' 
#' @export
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}