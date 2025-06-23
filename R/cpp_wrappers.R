#' @useDynLib manifoldhrf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Compute Pairwise Distances
#' 
#' Computes pairwise Euclidean distances between columns of a matrix.
#' 
#' @param M Matrix with features in rows and samples in columns
#' @return Distance matrix
#' @keywords internal
#' @export
pairwise_distances_cpp <- function(M) {
  .Call(`_manifoldhrf_pairwise_distances_cpp`, M)
}

#' K-Nearest Neighbors Search
#' 
#' Finds k nearest neighbors for query points in data.
#' 
#' @param data Data matrix (features × samples)
#' @param query Query matrix (features × samples)
#' @param k Number of nearest neighbors
#' @return List with idx (neighbor indices) and dist (distances)
#' @keywords internal
#' @export
knn_search_cpp <- function(data, query, k) {
  .Call(`_manifoldhrf_knn_search_cpp`, data, query, k)
}

