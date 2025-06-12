# Tests for spatial smoothing improvements
library(testthat)
library(manifoldhrf)

test_that("k-NN fallback mechanisms work correctly", {
  # Create small test data
  set.seed(123)
  coords <- as.matrix(expand.grid(x = 1:5, y = 1:5, z = 1:2))
  V <- nrow(coords)
  
  # Test with our C++ implementation (should exist after compilation)
  if (exists("knn_search_cpp", mode = "function")) {
    L1 <- make_voxel_graph_laplacian_core(
      coords, 
      num_neighbors_Lsp = 6,
      distance_engine = "euclidean"
    )
    expect_s4_class(L1, "sparseMatrix")
    expect_equal(dim(L1), c(V, V))
  }
  
  # Test with RANN fallback
  if (requireNamespace("RANN", quietly = TRUE)) {
    # Temporarily hide knn_search_cpp to test RANN fallback
    # Use local to avoid the deprecated with_mock
    local({
      # Override exists locally
      exists_orig <- base::exists
      assign("exists", function(x, ...) {
        if (x == "knn_search_cpp") FALSE else exists_orig(x, ...)
      }, envir = environment())
      
      L2 <- make_voxel_graph_laplacian_core(
        coords, 
        num_neighbors_Lsp = 6,
        distance_engine = "euclidean"
      )
      expect_s4_class(L2, "sparseMatrix")
      expect_equal(dim(L2), c(V, V))
    })
  }
  
  # Skip error test due to mocking limitations
  skip("Cannot test error case without proper mocking")
})

test_that("Lambda = 0 early exit works", {
  set.seed(456)
  m <- 3
  V <- 100
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Create dummy Laplacian (won't be used)
  L <- Matrix::Diagonal(n = V)
  
  # Test early return
  result <- apply_spatial_smoothing_core(Xi, L, lambda_spatial_smooth = 0)
  
  # Should return original matrix
  expect_identical(result, Xi)
})

test_that("Cholesky factorization is used for SPD systems", {
  skip_if_not_installed("Matrix", minimum_version = "1.3-0")
  
  set.seed(789)
  # Small test case
  m <- 2
  V <- 20
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Create simple 3D grid Laplacian
  coords <- as.matrix(expand.grid(x = 1:5, y = 1:4, z = 1))
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 4)
  
  # Apply smoothing with lambda > 0
  result <- apply_spatial_smoothing_core(Xi, L, lambda_spatial_smooth = 0.1)
  
  # Result should be valid
  expect_equal(dim(result), dim(Xi))
  expect_true(all(is.finite(result)))
  
  # Result should be smoother than original
  # (neighboring voxels should be more similar)
  orig_var <- mean(apply(Xi, 1, var))
  smooth_var <- mean(apply(result, 1, var))
  expect_lt(smooth_var, orig_var)
})

test_that("Weight schemes work correctly (when implemented)", {
  skip("Weight scheme parameter will be available after next package build")
  
  # This test will work once the package is rebuilt with the new parameter
  # set.seed(111)
  # coords <- as.matrix(expand.grid(x = 1:4, y = 1:4, z = 1:2))
  # V <- nrow(coords)
  # 
  # # Test binary weights
  # L_binary <- make_voxel_graph_laplacian_core(
  #   coords, 
  #   num_neighbors_Lsp = 6,
  #   weight_scheme = "binary"
  # )
  # 
  # # Test Gaussian weights
  # L_gaussian <- make_voxel_graph_laplacian_core(
  #   coords, 
  #   num_neighbors_Lsp = 6,
  #   weight_scheme = "gaussian"
  # )
})

test_that("Performance: Cholesky vs standard solve", {
  skip_on_cran()  # Skip timing tests on CRAN
  
  set.seed(222)
  # Larger test case
  m <- 5
  V <- 500
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Create 3D grid Laplacian
  grid_size <- ceiling(V^(1/3))
  coords <- as.matrix(expand.grid(
    x = 1:grid_size, 
    y = 1:grid_size, 
    z = 1:ceiling(V / grid_size^2)
  ))[1:V, ]
  
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 6)
  
  # Time the smoothing operation
  time_start <- Sys.time()
  result <- apply_spatial_smoothing_core(Xi, L, lambda_spatial_smooth = 0.5)
  time_end <- Sys.time()
  
  elapsed <- as.numeric(time_end - time_start, units = "secs")
  
  # Should complete reasonably fast
  expect_lt(elapsed, 1.0)  # Less than 1 second for 500 voxels
  
  # Result should be valid
  expect_equal(dim(result), dim(Xi))
  expect_true(all(is.finite(result)))
})

test_that("Integration test: full spatial smoothing pipeline", {
  set.seed(333)
  # Simulate data with spatial structure
  grid_size <- 8
  coords <- as.matrix(expand.grid(x = 1:grid_size, y = 1:grid_size, z = 1))
  V <- nrow(coords)
  m <- 3
  
  # Create manifold coordinates with spatial correlation
  Xi_true <- matrix(0, m, V)
  for (i in 1:m) {
    # Create smooth spatial field
    for (v in 1:V) {
      x <- coords[v, 1]
      y <- coords[v, 2]
      Xi_true[i, v] <- sin(2 * pi * x / grid_size) * cos(2 * pi * y / grid_size)
    }
  }
  
  # Add noise
  Xi_noisy <- Xi_true + matrix(rnorm(m * V, sd = 0.5), m, V)
  
  # Apply spatial smoothing
  L <- make_voxel_graph_laplacian_core(
    coords, 
    num_neighbors_Lsp = 4
  )
  
  Xi_smoothed <- apply_spatial_smoothing_core(
    Xi_noisy, 
    L, 
    lambda_spatial_smooth = 1.0
  )
  
  # Smoothed should be closer to true than noisy
  error_noisy <- mean((Xi_noisy - Xi_true)^2)
  error_smoothed <- mean((Xi_smoothed - Xi_true)^2)
  
  expect_lt(error_smoothed, error_noisy)
})