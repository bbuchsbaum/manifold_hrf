library(manifoldhrf)

compute_edge_weights_old <- function(Xi_matrix, voxel_coords, edge_threshold = 2) {
  V <- ncol(Xi_matrix)
  edge_weights <- rep(1, V)
  for (v in seq_len(V)) {
    voxel_v_matrix <- matrix(rep(voxel_coords[v, ], each = V), nrow = V, ncol = 3)
    distances <- sqrt(rowSums((voxel_coords - voxel_v_matrix)^2))
    neighbors <- which(distances > 0 & distances < 2)
    if (length(neighbors) > 0) {
      xi_v <- Xi_matrix[, v]
      xi_neighbors <- Xi_matrix[, neighbors, drop = FALSE]
      gradient_mag <- mean(colSums((xi_neighbors - xi_v)^2))
      if (gradient_mag > edge_threshold) {
        edge_weights[v] <- exp(-gradient_mag / edge_threshold)
      }
    }
  }
  edge_weights
}


test_that("optimized compute_edge_weights is faster", {
  skip_if_not_installed("microbenchmark")

  set.seed(123)
  coords <- as.matrix(expand.grid(x = 1:10, y = 1:10, z = 1:3))
  V <- nrow(coords)
  Xi <- matrix(rnorm(3 * V), 3, V)

  bm <- microbenchmark::microbenchmark(
    old = compute_edge_weights_old(Xi, coords),
    new = compute_edge_weights(Xi, coords),
    times = 3L
  )

  median_old <- median(bm$time[bm$expr == "old"])
  median_new <- median(bm$time[bm$expr == "new"])

  expect_lt(median_new, median_old)
})
