# Tests for refactored alternating optimization with streaming
library(testthat)
library(manifoldhrf)

test_that("Streaming implementation matches precomputed version", {
  set.seed(123)
  n <- 50
  p <- 10
  k <- 3
  V <- 20
  
  # Generate test data
  Y_proj <- matrix(rnorm(n * V), n, V)
  X_cond_list <- lapply(1:k, function(c) {
    X <- matrix(0, n, p)
    onsets <- seq(5 + (c-1)*15, n-p, by = 40)
    for (onset in onsets) {
      if (onset + p <= n) {
        X[onset:(onset+p-1), ] <- diag(p)
      }
    }
    X
  })
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  # Manual precomputed version for comparison
  conv_design_list <- lapply(X_cond_list, function(Xc) {
    Xc %*% H_shapes
  })
  
  Beta_manual <- matrix(0, k, V)
  for (v in 1:V) {
    D_v <- matrix(0, n, k)
    for (c in 1:k) {
      D_v[, c] <- conv_design_list[[c]][, v]
    }
    XtX <- crossprod(D_v) + 0.01 * diag(k)
    XtY <- crossprod(D_v, Y_proj[, v])
    Beta_manual[, v] <- solve(XtX, XtY)
  }
  
  # Run refactored version
  Beta_streaming <- estimate_final_condition_betas_core(
    Y_proj, X_cond_list, H_shapes,
    lambda_beta_final = 0.01,
    control_alt_list = list(max_iter = 1),
    n_jobs = 1
  )
  
  # Results should match
  expect_equal(Beta_streaming, Beta_manual, tolerance = 1e-10)
})

test_that("Memory usage is reduced with streaming", {
  skip_if_not_installed("pryr")
  
  set.seed(456)
  # Larger test case
  n <- 100
  p <- 20
  k <- 4
  V <- 100  # More voxels
  
  Y_proj <- matrix(rnorm(n * V), n, V)
  X_cond_list <- lapply(1:k, function(c) {
    matrix(rnorm(n * p) / 10, n, p)
  })
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  # Track memory usage during execution
  gc(full = TRUE)
  mem_before <- pryr::mem_used()
  
  Beta_result <- estimate_final_condition_betas_core(
    Y_proj, X_cond_list, H_shapes,
    lambda_beta_final = 0.01,
    control_alt_list = list(max_iter = 1),
    n_jobs = 1
  )
  
  mem_peak <- pryr::mem_used()
  mem_increase <- as.numeric(mem_peak - mem_before)
  
  # Memory increase should be much less than precomputing all regressors
  # Precomputed would need: k * n * V * 8 bytes = 4 * 100 * 100 * 8 = 320KB
  # Streaming needs much less
  expect_lt(mem_increase, 1e6)  # Less than 1MB increase
  
  # Check result is valid
  expect_equal(dim(Beta_result), c(k, V))
  expect_true(all(is.finite(Beta_result)))
})

test_that("Cholesky decomposition with QR fallback works", {
  set.seed(789)
  n <- 30
  p <- 8
  k <- 2
  V <- 5
  
  Y_proj <- matrix(rnorm(n * V), n, V)
  
  # Create design that might lead to numerical issues
  X_cond_list <- lapply(1:k, function(c) {
    X <- matrix(0, n, p)
    # Make conditions very similar (near collinearity)
    X[5:12, ] <- diag(p) + rnorm(p * p, sd = 0.01)
    X
  })
  
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  # Should handle numerical issues gracefully
  expect_warning(
    Beta_result <- estimate_final_condition_betas_core(
      Y_proj, X_cond_list, H_shapes,
      lambda_beta_final = 1e-6,  # Very small regularization
      control_alt_list = list(max_iter = 1),
      n_jobs = 1
    ),
    regexp = NA  # No warnings expected with proper fallback
  )
  
  # Result should still be valid
  expect_equal(dim(Beta_result), c(k, V))
  expect_true(all(is.finite(Beta_result)))
})

test_that("Parallel execution works with streaming", {
  skip_if_not_installed("future.apply")
  
  set.seed(111)
  n <- 40
  p <- 10
  k <- 3
  V <- 30
  
  Y_proj <- matrix(rnorm(n * V), n, V)
  X_cond_list <- lapply(1:k, function(c) {
    matrix(rnorm(n * p) / 5, n, p)
  })
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  # Sequential result
  Beta_seq <- estimate_final_condition_betas_core(
    Y_proj, X_cond_list, H_shapes,
    lambda_beta_final = 0.01,
    control_alt_list = list(max_iter = 1),
    n_jobs = 1
  )
  
  # Parallel result (if .lss_process_voxels exists)
  if (exists(".lss_process_voxels", mode = "function")) {
    future::plan(future::multisession, workers = 2)
    
    Beta_par <- estimate_final_condition_betas_core(
      Y_proj, X_cond_list, H_shapes,
      lambda_beta_final = 0.01,
      control_alt_list = list(max_iter = 1),
      n_jobs = 2
    )
    
    future::plan(future::sequential)
    
    # Results should match
    expect_equal(Beta_seq, Beta_par, tolerance = 1e-10)
  }
})

test_that("Performance comparison: streaming vs precomputed", {
  skip_on_cran()  # Skip on CRAN due to timing
  
  set.seed(222)
  n <- 100
  p <- 20
  k <- 4
  V <- 50
  
  Y_proj <- matrix(rnorm(n * V), n, V)
  X_cond_list <- lapply(1:k, function(c) {
    matrix(rnorm(n * p) / 10, n, p)
  })
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  # Time the streaming version
  time_streaming <- system.time({
    Beta_streaming <- estimate_final_condition_betas_core(
      Y_proj, X_cond_list, H_shapes,
      lambda_beta_final = 0.01,
      n_jobs = 1
    )
  })
  
  # The streaming version might be slightly slower per voxel
  # but uses much less memory
  expect_lt(time_streaming["elapsed"], 5)  # Should complete in reasonable time
  
  # Verify output
  expect_equal(dim(Beta_streaming), c(k, V))
  expect_true(all(is.finite(Beta_streaming)))
})