# Fallback distance utilities

#' Pairwise Euclidean distances
#'
#' If the Rcpp implementation is unavailable, this R version computes
#' pairwise Euclidean distances between the columns of a matrix.
#'
#' @param X Numeric matrix with observations in columns.
#' @return Matrix of Euclidean distances between columns of \code{X}.
#' @keywords internal
pairwise_distances_cpp <- function(X) {
  n_cols <- ncol(X)
  
  # Check memory requirements before computing distances
  tryCatch({
    check_distance_memory(n_cols)
  }, error = function(e) {
    stop("Cannot compute distance matrix: ", e$message, call. = FALSE)
  })
  
  as.matrix(dist(t(X)))
}
