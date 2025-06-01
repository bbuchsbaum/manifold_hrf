# Tests for Core Spatial Smoothing Functions (Component 2)
# Tests for MHRF-CORE-SPSMOOTH-01 and MHRF-CORE-SPSMOOTH-02

test_that("make_voxel_graph_laplacian_core works for simple grid", {
  # Create a simple 3x3x1 grid
  coords <- expand.grid(x = 1:3, y = 1:3, z = 1)
  voxel_coords <- as.matrix(coords)
  
  # Test with 4-neighbor connectivity (face-connected in 2D)
  L_sparse <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 4)
  
  # Check dimensions
  expect_equal(dim(L_sparse), c(9, 9))
  
  # Check that it's sparse
  expect_true(inherits(L_sparse, "sparseMatrix"))
  
  # Check Laplacian properties
  # Row sums should be zero (or very close due to numerical precision)
  row_sums <- Matrix::rowSums(L_sparse)
  expect_true(all(abs(row_sums) < 1e-10))
  
  # Diagonal should be positive (degrees)
  expect_true(all(Matrix::diag(L_sparse) > 0))
  
  # Off-diagonal should be non-positive
  L_dense <- as.matrix(L_sparse)
  diag(L_dense) <- 0
  expect_true(all(L_dense <= 0))
  
  # Check symmetry
  expect_true(Matrix::isSymmetric(L_sparse))
})

test_that("make_voxel_graph_laplacian_core handles different neighbor counts", {
  # Create a 5x5x5 cube
  coords <- expand.grid(x = 1:5, y = 1:5, z = 1:5)
  voxel_coords <- as.matrix(coords)
  
  # Test different neighbor counts
  L_6 <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 6)
  L_18 <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 18)
  L_26 <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 26)
  
  # More neighbors should lead to higher connectivity (more non-zero entries)
  expect_lt(Matrix::nnzero(L_6), Matrix::nnzero(L_18))
  expect_lt(Matrix::nnzero(L_18), Matrix::nnzero(L_26))
  
  # All should be valid Laplacians
  for (L in list(L_6, L_18, L_26)) {
    expect_true(all(abs(Matrix::rowSums(L)) < 1e-10))
    expect_true(Matrix::isSymmetric(L))
  }
})

test_that("make_voxel_graph_laplacian_core validates inputs", {
  # Valid coordinates
  coords <- matrix(1:15, 5, 3)
  
  # Test non-matrix input
  expect_error(
    make_voxel_graph_laplacian_core(data.frame(coords), 6),
    "must be a matrix"
  )
  
  # Test wrong number of columns
  expect_error(
    make_voxel_graph_laplacian_core(matrix(1:10, 5, 2), 6),
    "must have exactly 3 columns"
  )
  
  # Test invalid neighbor count
  expect_error(
    make_voxel_graph_laplacian_core(coords, -1),
    "must be a positive integer"
  )
  
  expect_error(
    make_voxel_graph_laplacian_core(coords, 2.5),
    "must be a positive integer"
  )
  
  # Test too few voxels
  expect_error(
    make_voxel_graph_laplacian_core(matrix(1:3, 1, 3), 6),
    "must have at least 2 rows"
  )
})

test_that("make_voxel_graph_laplacian_core handles edge cases", {
  # Test with just 2 voxels
  coords <- matrix(c(1, 2, 1, 1, 1, 1), 2, 3)
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 1)
  
  # Should be a 2x2 matrix with specific structure
  expect_equal(dim(L), c(2, 2))
  expect_equal(as.matrix(L), matrix(c(1, -1, -1, 1), 2, 2))
  
  # Test requesting more neighbors than available
  coords <- matrix(rnorm(15), 5, 3)
  expect_warning(
    L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 10),
    "only 4 other voxels available"
  )
  
  # Should still produce valid Laplacian
  expect_true(all(abs(Matrix::rowSums(L)) < 1e-10))
})

test_that("ann_euclidean falls back when RcppHNSW absent", {
  coords <- matrix(runif(30), 10, 3)

  if (!requireNamespace("RcppHNSW", quietly = TRUE)) {
    expect_warning(
      L_ann <- make_voxel_graph_laplacian_core(
        coords, num_neighbors_Lsp = 3, distance_engine = "ann_euclidean"
      ),
      "falling back"
    )
    L_exact <- make_voxel_graph_laplacian_core(
      coords, num_neighbors_Lsp = 3, distance_engine = "euclidean"
    )
    expect_equal(L_ann, L_exact)
  } else {
    L_ann <- make_voxel_graph_laplacian_core(
      coords, num_neighbors_Lsp = 3, distance_engine = "ann_euclidean"
    )
    expect_true(inherits(L_ann, "Matrix"))
  }
})

test_that("apply_spatial_smoothing_core works correctly", {
  # Create test data
  set.seed(123)
  m <- 4  # manifold dimensions
  V <- 25 # voxels (5x5 grid)
  
  # Create grid coordinates
  coords <- expand.grid(x = 1:5, y = 1:5, z = 1)
  voxel_coords <- as.matrix(coords)
  
  # Create Laplacian
  L_sparse <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 4)
  
  # Create manifold coordinates with spatial structure
  Xi_ident <- matrix(0, m, V)
  for (j in seq_len(m)) {
    # Add smooth spatial pattern
    Xi_ident[j, ] <- sin(coords$x / 2) + cos(coords$y / 2)
    # Add noise
    Xi_ident[j, ] <- Xi_ident[j, ] + rnorm(V, sd = 0.5)
  }
  
  # Apply smoothing
  Xi_smoothed <- apply_spatial_smoothing_core(Xi_ident, L_sparse, lambda_spatial_smooth = 0.5)

  # Check dimensions
  expect_equal(dim(Xi_smoothed), c(m, V))
  expect_true(is.matrix(Xi_smoothed))
  
  # Smoothed should have less variance than original
  for (j in seq_len(m)) {
    expect_lt(var(Xi_smoothed[j, ]), var(Xi_ident[j, ]))
  }
  
  # Smoothed coordinates should be spatially smoother
  # Check by computing Laplacian quadratic form: x'Lx
  for (j in seq_len(m)) {
    roughness_original <- as.numeric(t(Xi_ident[j, ]) %*% L_sparse %*% Xi_ident[j, ])
    roughness_smoothed <- as.numeric(t(Xi_smoothed[j, ]) %*% L_sparse %*% Xi_smoothed[j, ])
    expect_lt(roughness_smoothed, roughness_original)
  }
})

test_that("apply_spatial_smoothing_core with lambda=0 returns original", {
  # With lambda=0, should return original coordinates
  set.seed(456)
  m <- 3
  V <- 20
  
  Xi_ident <- matrix(rnorm(m * V), m, V)
  
  coords <- matrix(rnorm(V * 3), V, 3)
  L_sparse <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 5)
  
  # Apply with lambda = 0
  Xi_smoothed <- apply_spatial_smoothing_core(Xi_ident, L_sparse, lambda_spatial_smooth = 0)
  
  # Should be identical to original
  expect_equal(Xi_smoothed, Xi_ident, tolerance = 1e-10)
})

test_that("apply_spatial_smoothing_core matches column-wise solve", {
  set.seed(42)
  m <- 3
  V <- 15
  Xi_ident <- matrix(rnorm(m * V), m, V)
  coords <- matrix(rnorm(V * 3), V, 3)
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 5)
  lambda <- 0.2

  Xi_new <- apply_spatial_smoothing_core(Xi_ident, L, lambda)

  A <- Matrix::Diagonal(V) + lambda * L
  Xi_expected <- matrix(0, m, V)
  for (j in 1:m) {
    Xi_expected[j, ] <- as.vector(Matrix::solve(A, Xi_ident[j, ]))
  }

  expect_equal(Xi_new, Xi_expected)
  expect_true(is.matrix(Xi_new))
})

test_that("apply_spatial_smoothing_core validates inputs", {
  # Valid inputs
  Xi <- matrix(1:20, 4, 5)
  coords <- matrix(rnorm(15), 5, 3)
  L <- make_voxel_graph_laplacian_core(coords, 3)
  
  # Test non-matrix Xi
  expect_error(
    apply_spatial_smoothing_core(data.frame(Xi), L, 0.1),
    "Xi_ident_matrix must be a matrix"
  )
  
  # Test non-sparse L
  expect_error(
    apply_spatial_smoothing_core(Xi, as.matrix(L), 0.1),
    "must be a sparse matrix"
  )
  
  # Test invalid lambda
  expect_error(
    apply_spatial_smoothing_core(Xi, L, -0.1),
    "must be a non-negative scalar"
  )
  
  expect_error(
    apply_spatial_smoothing_core(Xi, L, c(0.1, 0.2)),
    "must be a non-negative scalar"
  )
  
  # Test dimension mismatch
  Xi_bad <- matrix(1:24, 4, 6)  # Wrong number of voxels
  expect_error(
    apply_spatial_smoothing_core(Xi_bad, L, 0.1),
    "must be 6 x 6 to match Xi_ident_matrix"
  )
})

test_that("apply_spatial_smoothing_core handles zero dimensions", {
  # Case 1: zero manifold dimensions
  V <- 5
  Xi_zero_m <- matrix(numeric(0), nrow = 0, ncol = V)
  coords <- matrix(rnorm(V * 3), V, 3)
  L_sparse <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 3)
  result_m0 <- apply_spatial_smoothing_core(Xi_zero_m, L_sparse, lambda_spatial_smooth = 0.1)
  expect_equal(result_m0, Xi_zero_m)

  # Case 2: zero voxels
  Xi_zero_v <- matrix(numeric(0), nrow = 3, ncol = 0)
  L_empty <- Matrix::Matrix(0, 0, 0, sparse = TRUE)
  result_v0 <- apply_spatial_smoothing_core(Xi_zero_v, L_empty, lambda_spatial_smooth = 0.1)
  expect_equal(result_v0, Xi_zero_v)
})

test_that("apply_spatial_smoothing_core preserves mean", {
  # Smoothing should approximately preserve the mean of each manifold dimension
  set.seed(789)
  m <- 5
  V <- 50
  
  Xi_ident <- matrix(rnorm(m * V, mean = 2), m, V)
  
  # Create random coordinates
  coords <- matrix(runif(V * 3, 0, 10), V, 3)
  L_sparse <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 10)
  
  # Apply moderate smoothing
  Xi_smoothed <- apply_spatial_smoothing_core(Xi_ident, L_sparse, lambda_spatial_smooth = 0.2)
  
  # Check that means are preserved (approximately)
  for (j in seq_len(m)) {
    mean_original <- mean(Xi_ident[j, ])
    mean_smoothed <- mean(Xi_smoothed[j, ])
    expect_equal(mean_smoothed, mean_original, tolerance = 0.01)
  }
})

test_that("spatial smoothing integration test", {
  # Test the full spatial smoothing pipeline
  set.seed(321)
  
  # Create a 10x10x3 brain-like grid
  coords <- expand.grid(x = 1:10, y = 1:10, z = 1:3)
  voxel_coords <- as.matrix(coords)
  V <- nrow(voxel_coords)
  m <- 4
  
  # Create manifold coordinates with spatial clusters
  Xi_ident <- matrix(0, m, V)
  
  # Add different patterns to each manifold dimension
  Xi_ident[1, ] <- ifelse(coords$x < 5, 1, -1) + rnorm(V, sd = 0.3)
  Xi_ident[2, ] <- ifelse(coords$y < 5, 1, -1) + rnorm(V, sd = 0.3)
  Xi_ident[3, ] <- sin(coords$x / 2) * cos(coords$y / 2) + rnorm(V, sd = 0.3)
  Xi_ident[4, ] <- coords$z / 3 + rnorm(V, sd = 0.3)
  
  # Create Laplacian with 26-connectivity
  L_sparse <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 26)
  
  # Test different smoothing levels
  lambdas <- c(0, 0.1, 0.5, 1.0, 5.0)
  roughness_by_lambda <- matrix(0, m, length(lambdas))
  
  for (i in seq_along(lambdas)) {
    Xi_smoothed <- apply_spatial_smoothing_core(Xi_ident, L_sparse, lambdas[i])
    
    # Compute roughness for each dimension
    for (j in seq_len(m)) {
      roughness_by_lambda[j, i] <- as.numeric(
        t(Xi_smoothed[j, ]) %*% L_sparse %*% Xi_smoothed[j, ]
      )
    }
  }
  
  # Roughness should decrease with increasing lambda
  for (j in seq_len(m)) {
    expect_true(all(diff(roughness_by_lambda[j, ]) <= 0))
  }
  
  # Very high smoothing should make coordinates nearly constant
  Xi_very_smooth <- apply_spatial_smoothing_core(Xi_ident, L_sparse, lambda_spatial_smooth = 100)
  for (j in seq_len(m)) {
    # Variance should be very small
    expect_lt(var(Xi_very_smooth[j, ]), 0.01)
  }
})

test_that("symmetrization step retains sparse class", {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    M <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 3), x = 1, dims = c(3, 3))
    M_sym <- Matrix::pmax(M, Matrix::t(M))
    expect_true(inherits(M_sym, "sparseMatrix"))
  }
})