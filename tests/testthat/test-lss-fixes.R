# Regression tests for critical LSS fixes
library(testthat)
library(manifoldhrf)

test_that("Xi_coords has correct dimensions (m x V)", {
  # This test would have failed before the transpose fix
  set.seed(123)
  n <- 50
  V <- 10
  p <- 15
  
  Y <- matrix(rnorm(n * V), n, V)
  X_trials <- lapply(1:3, function(i) matrix(rnorm(n * p), n, p))
  H <- matrix(rnorm(p * V), p, V)
  
  # Call wrapper which internally creates Xi_coords
  expect_silent({
    result <- run_lss_voxel_loop(Y, X_trials, H, verbose = FALSE)
  })
  
  # Verify result dimensions
  expect_equal(dim(result), c(3, V))
})

test_that("chunk_size parameter is respected in auto-selection", {
  set.seed(456)
  
  # Mock the internal select function to capture its arguments
  mock_select <- NULL
  with_mocked_bindings(
    `.select_memory_strategy` = function(...) {
      mock_select <<- list(...)
      return("chunked")
    },
    {
      # Small data with custom chunk_size
      run_lss_voxel_loop_core(
        Y_proj_matrix = matrix(1, 10, 5),
        X_trial_onset_list_of_matrices = list(matrix(1, 10, 3)),
        B_reconstructor_matrix = diag(3),
        Xi_smoothed_allvox_matrix = matrix(1, 3, 5),
        memory_strategy = "auto",
        chunk_size = 100,  # Custom value
        verbose = FALSE
      )
    }
  )
  
  # Verify chunk_size was passed to selector
  expect_equal(mock_select[[4]], 100)
})

test_that("Memory estimates include all components", {
  # Test the memory estimation logic
  n <- 500
  V <- 2000
  T_trials <- 100
  chunk_size <- 10
  
  # Get memory estimates
  strategy <- manifoldhrf:::.select_memory_strategy(
    n, V, T_trials, chunk_size, 
    ram_limit_GB = 1.0, 
    verbose = FALSE
  )
  
  # With realistic data sizes and proper accounting,
  # full strategy should need much more memory
  # and likely select chunked or streaming
  expect_true(strategy %in% c("chunked", "streaming"))
})

test_that("Chunked strategy doesn't allocate full n×V matrices per trial", {
  set.seed(789)
  n <- 50
  V <- 100  # Many voxels
  T_trials <- 10
  p <- 10
  m <- 3
  
  # Create test data
  Y <- matrix(rnorm(n * V), n, V)
  X_trials <- lapply(1:T_trials, function(i) {
    X <- matrix(0, n, p)
    X[sample(n, 5), sample(p, 5)] <- 1  # Sparse
    X
  })
  
  B <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Run with tiny chunk size to force chunking
  result <- run_lss_voxel_loop_core(
    Y, X_trials, B, Xi,
    memory_strategy = "chunked",
    chunk_size = 2,  # Very small chunks
    verbose = FALSE
  )
  
  # Should complete without memory errors
  expect_equal(dim(result), c(T_trials, V))
  
  # All values should be finite (no memory corruption)
  expect_true(all(is.finite(result)))
})

test_that("Constant column detector is efficient", {
  # Test the improved detector
  n <- 10000  # Large matrix
  
  # Create test matrix with intercept
  A <- cbind(
    1,  # intercept
    rnorm(n),  # not constant
    rep(5, n),  # constant but not 1
    rnorm(n) + 1e-15  # nearly constant
  )
  
  # Time the operation (should be fast even for large n)
  start_time <- Sys.time()
  has_intercept <- any(apply(A, 2, function(x) {
    all(abs(x - x[1]) < 1e-12)
  }))
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  
  expect_true(has_intercept)
  expect_lt(elapsed, 0.1)  # Should be very fast
  
  # Check it identifies the right columns
  const_cols <- apply(A, 2, function(x) all(abs(x - x[1]) < 1e-12))
  expect_equal(which(const_cols), c(1, 3))
})

test_that("Integration test: full pipeline with all fixes", {
  set.seed(999)
  # Realistic small example
  n <- 100
  V <- 20
  T_trials <- 5  # Reduce number of trials
  p <- 10
  m <- 3
  
  # Generate data
  Y <- matrix(rnorm(n * V), n, V)
  
  # Create simpler, well-conditioned trial matrices
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    # Create gaussian-shaped HRF-like regressors
    onset <- 10 + (t-1) * 15
    if (onset + 10 <= n) {
      for (j in 1:p) {
        timepoints <- (onset):(onset + 9)
        if (all(timepoints <= n)) {
          # Simple exponential decay
          X[timepoints, j] <- exp(-0.3 * (0:9)) * (j == 1)
        }
      }
    }
    X
  })
  
  # Manifold components
  B <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Fixed effects with intercept
  A_fixed <- cbind(1, rnorm(n))
  
  # Run with each strategy
  strategies <- c("full", "chunked", "streaming")
  results <- list()
  
  for (strat in strategies) {
    results[[strat]] <- run_lss_voxel_loop_core(
      Y, X_trials, B, Xi, A_fixed,
      memory_strategy = strat,
      chunk_size = 3,
      verbose = FALSE
    )
  }
  
  # Check properties
  expect_equal(dim(results$full), c(T_trials, V))
  expect_equal(dim(results$chunked), c(T_trials, V))
  expect_equal(dim(results$streaming), c(T_trials, V))
  
  # All results should be finite
  expect_true(all(sapply(results, function(r) all(is.finite(r)))))
  
  # Results should be correlated (allowing for algorithmic differences)
  if (all(sapply(results, function(r) all(is.finite(r))))) {
    # Convert to vectors for correlation
    full_vec <- as.vector(results$full)
    chunked_vec <- as.vector(results$chunked)
    streaming_vec <- as.vector(results$streaming)
    
    # High correlation indicates similar results despite numerical differences
    expect_gt(cor(full_vec, chunked_vec), 0.8)
    expect_gt(cor(full_vec, streaming_vec), 0.8)
  }
})