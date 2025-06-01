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
  as.matrix(dist(t(X)))
}
