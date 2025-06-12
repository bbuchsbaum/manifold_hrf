# Test parallel execution in LSS
library(testthat)
library(manifoldhrf)

test_that("Parallel execution produces identical results to sequential", {
  # Skip if future not available
  skip_if_not_installed("future.apply")
  
  set.seed(123)
  n <- 50
  V <- 20  # Small enough to be fast
  T_trials <- 5
  p <- 10
  m <- 3
  
  # Create test data
  Y <- matrix(rnorm(n * V), n, V)
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * 9
    if (onset + p <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })
  
  B <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Run sequential (n_cores = 1)
  result_seq <- run_lss_voxel_loop_core(
    Y, X_trials, B, Xi,
    memory_strategy = "full",
    n_cores = 1,
    progress = FALSE,
    verbose = FALSE
  )
  
  # Set up parallel backend
  future::plan(future::multisession, workers = 2)
  
  # Run parallel (n_cores = 2)
  result_par <- run_lss_voxel_loop_core(
    Y, X_trials, B, Xi,
    memory_strategy = "full",
    n_cores = 2,
    progress = FALSE,
    verbose = FALSE
  )
  
  # Reset to sequential
  future::plan(future::sequential)
  
  # Results should be identical
  expect_equal(result_seq, result_par, tolerance = 1e-10)
})

test_that("Parallel execution works for all strategies", {
  skip_if_not_installed("future.apply")
  
  set.seed(456)
  n <- 40
  V <- 15
  T_trials <- 6
  p <- 8
  m <- 3
  
  Y <- matrix(rnorm(n * V), n, V)
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    X[sample(n, 5), sample(p, 5)] <- 1
    X
  })
  
  B <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
  Xi <- matrix(rnorm(m * V), m, V)
  
  strategies <- c("full", "chunked", "streaming")
  
  # Set up parallel
  future::plan(future::multisession, workers = 2)
  
  for (strat in strategies) {
    result <- run_lss_voxel_loop_core(
      Y, X_trials, B, Xi,
      memory_strategy = strat,
      chunk_size = 3,
      n_cores = 2,
      progress = FALSE,
      verbose = FALSE
    )
    
    expect_equal(dim(result), c(T_trials, V))
    expect_true(all(is.finite(result)))
  }
  
  future::plan(future::sequential)
})

test_that("Progress reporting works in parallel mode", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("progressr")
  
  set.seed(789)
  n <- 30
  V <- 10
  T_trials <- 3
  p <- 5
  m <- 2
  
  Y <- matrix(rnorm(n * V), n, V)
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    X[5:(5+p-1), ] <- diag(p)
    X
  })
  
  B <- diag(p)  # Should match trial matrix columns (p)
  Xi <- matrix(rnorm(p * V), p, V)  # Should match B columns
  
  # Set up progress handling
  progressr::with_progress({
    result <- run_lss_voxel_loop_core(
      Y, X_trials, B, Xi,
      memory_strategy = "streaming",
      n_cores = 2,
      progress = TRUE,
      verbose = FALSE
    )
  })
  
  expect_equal(dim(result), c(T_trials, V))
})

test_that("Parallel execution handles edge cases", {
  skip_if_not_installed("future.apply")
  
  set.seed(999)
  
  # Edge case: more cores than voxels
  Y_small <- matrix(rnorm(20 * 3), 20, 3)  # Only 3 voxels
  X_trials <- list(diag(5), diag(5))
  B <- matrix(1, 5, 1)
  Xi_small <- matrix(1, 1, 3)
  
  future::plan(future::multisession, workers = 4)
  
  result <- run_lss_voxel_loop_core(
    Y_small, X_trials, B, Xi_small,
    memory_strategy = "full",
    n_cores = 4,  # More cores than voxels
    progress = FALSE,
    verbose = FALSE
  )
  
  expect_equal(dim(result), c(2, 3))
  
  future::plan(future::sequential)
})

test_that("Worker function isolation is maintained", {
  skip_if_not_installed("future.apply")
  
  # Test that parallel workers don't interfere with each other
  set.seed(111)
  n <- 40
  V <- 50
  T_trials <- 4
  p <- 10
  m <- 3
  
  # Create data with known pattern
  Y <- matrix(0, n, V)
  for (v in 1:V) {
    Y[, v] <- sin(2 * pi * (1:n) / n + v/10)  # Different phase per voxel
  }
  
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    start_row <- (t-1) * 8 + 1  # More conservative spacing
    end_row <- min(start_row + p - 1, n)  # Don't exceed matrix bounds
    if (end_row >= start_row) {
      rows_available <- end_row - start_row + 1
      X[start_row:end_row, 1:rows_available] <- diag(rows_available)
    }
    X
  })
  
  B <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Run with multiple cores
  future::plan(future::multisession, workers = 3)
  
  result_parallel <- run_lss_voxel_loop_core(
    Y, X_trials, B, Xi,
    memory_strategy = "streaming",
    n_cores = 3,
    progress = FALSE,
    verbose = FALSE
  )
  
  future::plan(future::sequential)
  
  # Each voxel should have been processed independently
  # Check that results vary across voxels (not all identical)
  voxel_vars <- apply(result_parallel, 2, var)
  expect_true(var(voxel_vars) > 0)  # Voxels should have different patterns
})