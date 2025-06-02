# Tests for Core Alternating Optimization Functions (Component 4)
# Tests for MHRF-CORE-ALTOPT-01

test_that("estimate_final_condition_betas_core works correctly", {
  # Create test data
  set.seed(123)
  n <- 200  # timepoints
  p <- 25   # HRF length
  k <- 3    # conditions
  V <- 50   # voxels
  
  # Create condition design matrices with different event patterns
  X_cond_list <- list()
  for (c in 1:k) {
    X <- matrix(0, n, p)
    # Different onset patterns for each condition
    if (c == 1) {
      onsets <- seq(10, n-p, by = 40)  # Regular
    } else if (c == 2) {
      onsets <- seq(20, n-p, by = 50)  # Offset regular
    } else {
      onsets <- sort(sample(30:(n-p-10), 5))  # Random
    }
    
    for (onset in onsets) {
      # Simple FIR design
      for (j in 1:p) {
        if (onset + j - 1 <= n) {
          X[onset + j - 1, j] <- 1
        }
      }
    }
    X_cond_list[[c]] <- X
  }
  
  # Create HRF shapes with some variation
  H_shapes <- matrix(0, p, V)
  t_hrf <- seq(0, p-1) * 0.5
  for (v in 1:V) {
    # Gamma-like HRF with variation
    peak_time <- 5 + rnorm(1, 0, 0.5)
    h <- dgamma(t_hrf - 2, shape = peak_time, rate = 1)
    H_shapes[, v] <- h / sum(h)
  }
  
  # Create data with known condition betas
  true_betas <- matrix(0, k, V)
  true_betas[1, ] <- rnorm(V, mean = 1.5, sd = 0.3)  # Condition 1 positive
  true_betas[2, ] <- rnorm(V, mean = -0.5, sd = 0.3) # Condition 2 negative
  true_betas[3, ] <- rnorm(V, mean = 0, sd = 0.5)    # Condition 3 mixed
  
  # Generate Y data
  Y_proj <- matrix(0, n, V)
  for (v in 1:V) {
    signal <- rep(0, n)
    for (c in 1:k) {
      regressor <- X_cond_list[[c]] %*% H_shapes[, v]
      signal <- signal + true_betas[c, v] * regressor
    }
    Y_proj[, v] <- signal + rnorm(n, sd = 0.2)
  }
  
  # Estimate final betas
  Beta_estimated <- estimate_final_condition_betas_core(
    Y_proj, X_cond_list, H_shapes,
    lambda_beta_final = 0.01,
    control_alt_list = list(max_iter = 1)
  )
  
  # Check output
  expect_equal(dim(Beta_estimated), c(k, V))
  
  # Check correlation with true betas
  cor_overall <- cor(as.vector(true_betas), as.vector(Beta_estimated))
  expect_gt(cor_overall, 0.7)
  
  # Check condition-wise correlations
  # At least 2 out of 3 conditions should show good correlation
  cor_conditions <- sapply(1:k, function(c) {
    cor(true_betas[c, ], Beta_estimated[c, ])
  })
  expect_gt(sum(cor_conditions > 0.5), 1)  # At least 2 conditions > 0.5
  expect_gt(mean(cor_conditions), 0.4)      # Average correlation reasonable
})

test_that("estimate_final_condition_betas_core validates inputs", {
  # Valid base inputs
  n <- 100
  p <- 20
  k <- 2
  V <- 30
  
  Y <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  H <- matrix(rnorm(p * V), p, V)
  
  # Test non-matrix Y
  expect_error(
    estimate_final_condition_betas_core(data.frame(Y), X_list, H),
    "Y_proj_matrix must be a matrix"
  )
  
  # Test non-list X
  expect_error(
    estimate_final_condition_betas_core(Y, X_list[[1]], H),
    "X_condition_list_proj_matrices must be a list"
  )
  
  # Test empty condition list
  expect_error(
    estimate_final_condition_betas_core(Y, list(), H),
    "must contain at least one condition"
  )
  
  # Test non-matrix H
  expect_error(
    estimate_final_condition_betas_core(Y, X_list, data.frame(H)),
    "H_shapes_allvox_matrix must be a matrix"
  )
  
  # Test dimension mismatches
  H_bad <- matrix(rnorm(p * (V-1)), p, V-1)
  expect_error(
    estimate_final_condition_betas_core(Y, X_list, H_bad),
    "must have V columns"
  )
  
  # Test invalid lambda
  expect_error(
    estimate_final_condition_betas_core(Y, X_list, H, lambda_beta_final = -0.1),
    "must be a non-negative scalar"
  )
  
  # Test invalid control parameters
  expect_error(
    estimate_final_condition_betas_core(Y, X_list, H, 
                                      control_alt_list = "not a list"),
    "control_alt_list must be a list"
  )
  
  expect_error(
    estimate_final_condition_betas_core(Y, X_list, H,
                                      control_alt_list = list(max_iter = 0)),
    "max_iter must be a positive integer"
  )
  
  expect_error(
    estimate_final_condition_betas_core(Y, X_list, H,
                                      control_alt_list = list(rel_change_tol = -1)),
    "rel_change_tol must be a positive scalar"
  )
})

test_that("estimate_final_condition_betas_core handles edge cases", {
  set.seed(456)
  n <- 100
  p <- 20
  k <- 2
  V <- 10
  
  # Test with zero HRFs
  Y <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  H_zero <- matrix(0, p, V)
  
  # Should handle gracefully
  expect_warning(
    Beta_zero <- estimate_final_condition_betas_core(Y, X_list, H_zero),
    NA  # No specific warning expected, but function should handle it
  )
  
  # All betas should be zero
  expect_equal(sum(abs(Beta_zero)), 0)
  
  # Test with very small lambda (near zero)
  H_normal <- matrix(rnorm(p * V), p, V)
  Beta_small_lambda <- estimate_final_condition_betas_core(
    Y, X_list, H_normal,
    lambda_beta_final = 1e-10
  )
  
  # Should still produce valid output
  expect_false(any(is.na(Beta_small_lambda)))
  expect_false(any(is.infinite(Beta_small_lambda)))
})

test_that("estimate_final_condition_betas_core works with multiple iterations", {
  # Test the iterative refinement capability
  set.seed(789)
  n <- 150
  p <- 20
  k <- 2
  V <- 25
  
  # Create simple test case
  Y <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:k, function(i) {
    X <- matrix(0, n, p)
    # Add some events
    onsets <- seq(10, n-p, by = 30)
    for (onset in onsets) {
      X[onset:(onset+p-1), ] <- diag(p)
    }
    X
  })
  H <- matrix(abs(rnorm(p * V)), p, V)  # Positive HRFs
  
  # Run with multiple iterations
  Beta_iter1 <- estimate_final_condition_betas_core(
    Y, X_list, H,
    lambda_beta_final = 0.1,
    control_alt_list = list(max_iter = 1)
  )
  
  Beta_iter3 <- estimate_final_condition_betas_core(
    Y, X_list, H,
    lambda_beta_final = 0.1,
    control_alt_list = list(max_iter = 3)
  )
  
  # Both should produce valid results
  expect_equal(dim(Beta_iter1), c(k, V))
  expect_equal(dim(Beta_iter3), c(k, V))
  
  # With fixed HRFs, results should be identical
  # (since we're not actually updating HRFs between iterations)
  expect_equal(Beta_iter1, Beta_iter3, tolerance = 1e-8)
})

test_that("estimate_final_condition_betas_core recovers known signal patterns", {
  # Test with clearly separable conditions
  set.seed(321)
  n <- 300
  p <- 30
  k <- 3
  V <- 20
  
  # Create non-overlapping conditions
  X_cond_list <- list()
  
  # Condition 1: Early trials
  X_cond_list[[1]] <- matrix(0, n, p)
  for (onset in seq(10, 80, by = 20)) {
    X_cond_list[[1]][onset:(onset+p-1), ] <- diag(p)
  }
  
  # Condition 2: Middle trials
  X_cond_list[[2]] <- matrix(0, n, p)
  for (onset in seq(110, 180, by = 20)) {
    X_cond_list[[2]][onset:(onset+p-1), ] <- diag(p)
  }
  
  # Condition 3: Late trials
  X_cond_list[[3]] <- matrix(0, n, p)
  for (onset in seq(210, 280, by = 20)) {
    X_cond_list[[3]][onset:(onset+p-1), ] <- diag(p)
  }
  
  # Simple canonical HRF
  t_hrf <- seq(0, p-1) * 0.5
  h_canonical <- dgamma(t_hrf, shape = 6, rate = 1)
  h_canonical <- h_canonical / max(h_canonical)
  H_shapes <- matrix(rep(h_canonical, V), p, V)
  
  # Add slight variation
  H_shapes <- H_shapes + matrix(rnorm(p * V, sd = 0.05), p, V)
  
  # Known betas with clear patterns
  true_betas <- matrix(0, k, V)
  true_betas[1, 1:7] <- 2      # Condition 1 active in first third
  true_betas[2, 8:14] <- 1.5   # Condition 2 active in middle third
  true_betas[3, 15:20] <- 1    # Condition 3 active in last third
  
  # Generate clean data
  Y_proj <- matrix(0, n, V)
  for (v in 1:V) {
    for (c in 1:k) {
      regressor <- X_cond_list[[c]] %*% H_shapes[, v]
      Y_proj[, v] <- Y_proj[, v] + true_betas[c, v] * regressor
    }
    Y_proj[, v] <- Y_proj[, v] + rnorm(n, sd = 0.1)  # Small noise
  }
  
  # Estimate betas
  Beta_estimated <- estimate_final_condition_betas_core(
    Y_proj, X_cond_list, H_shapes,
    lambda_beta_final = 0.001
  )
  
  # Check activation patterns are recovered
  # Condition 1 should be highest in voxels 1-7
  for (v in 1:7) {
    expect_equal(which.max(Beta_estimated[, v]), 1)
  }
  
  # Condition 2 should be highest in voxels 8-14
  for (v in 8:14) {
    expect_equal(which.max(Beta_estimated[, v]), 2)
  }
  
  # Condition 3 should be highest in voxels 15-20
  for (v in 15:20) {
    expect_equal(which.max(Beta_estimated[, v]), 3)
  }
  
  # Correlation should be very high for this clean case
  cor_clean <- cor(as.vector(true_betas), as.vector(Beta_estimated))
  expect_gt(cor_clean, 0.9)
})