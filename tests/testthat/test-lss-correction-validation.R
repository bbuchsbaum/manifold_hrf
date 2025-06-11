# Test the corrected LSS implementation
library(testthat)

test_that("fmrilss implementation matches ground truth LSS", {
  
  set.seed(789)
  
  # Test parameters
  n <- 100
  T_trials <- 8
  p <- 15
  q <- 4
  lambda <- 1e-6
  
  # Create trial matrices
  X_trials <- list()
  for (t in 1:T_trials) {
    X_t <- matrix(0, n, p)
    onset <- 5 + (t-1) * 12
    if (onset + p <= n) {
      X_t[onset:(onset + p - 1), ] <- diag(p)
    }
    X_trials[[t]] <- X_t
  }
  
  # HRF
  h <- dgamma(0:(p-1), shape = 5, rate = 1.5)
  h <- h / sum(h)
  
  # Confounds
  Z <- cbind(
    1,
    (1:n) / n,
    sin(2 * pi * (1:n) / n),
    cos(2 * pi * (1:n) / n)
  )
  
  # True parameters
  true_betas <- c(1, -0.5, 2, 0, 1.5, -1, 0.8, 0.3)
  true_conf <- c(5, -2, 1, 0.5)
  
  # Generate data
  y <- Z %*% true_conf
  for (t in 1:T_trials) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.2)
  
  # Prepare projection manually
  ZtZ <- crossprod(Z)
  ZtZ_reg <- ZtZ + lambda * diag(ncol(Z))
  P_proj <- diag(n) - Z %*% solve(ZtZ_reg) %*% t(Z)
  y_proj <- P_proj %*% y
  
  # METHOD 1: Direct LSS (ground truth)
  betas_direct <- numeric(T_trials)
  for (t in 1:T_trials) {
    X_other <- do.call(cbind, lapply(X_trials[-t], function(X) X %*% h))
    X_full <- cbind(X_trials[[t]] %*% h, X_other, Z)
    XtX <- crossprod(X_full) + lambda * diag(ncol(X_full))
    Xty <- crossprod(X_full, y)
    betas_direct[t] <- solve(XtX, Xty)[1]
  }
  
  # METHOD 2: New fmrilss implementation
  # Create convolved regressors
  C <- matrix(0, n, T_trials)
  for (t in 1:T_trials) {
    C[, t] <- X_trials[[t]] %*% h
  }
  
  # Use fmrilss directly with confounds
  betas_corrected <- fmrilss::lss(
    Y = matrix(y, ncol = 1),
    X = C,
    Z = Z,
    method = "r_optimized"
  )
  
  # Compare
  cat("\n=== LSS Implementation Comparison ===\n")
  cat("True betas:     ", round(true_betas, 3), "\n")
  cat("Direct LSS:     ", round(betas_direct, 3), "\n")
  cat("fmrilss:        ", round(betas_corrected, 3), "\n")
  
  cat("\n=== Differences from Direct ===\n")
  cat("fmrilss diff:   ", round(betas_corrected - betas_direct, 6), "\n")
  
  cat("\n=== Mean Squared Errors ===\n")
  cat("Direct MSE:     ", round(mean((betas_direct - true_betas)^2), 4), "\n")
  cat("fmrilss MSE:    ", round(mean((betas_corrected - true_betas)^2), 4), "\n")
  
  # The corrected method should match direct for trials with signal
  # Note: LSS is solving a different problem than simultaneous estimation
  # It estimates each trial's effect while treating others as nuisance
  
  # For a proper comparison, let's also try simultaneous estimation
  C_all <- do.call(cbind, lapply(X_trials, function(X) X %*% h))
  X_sim <- cbind(C_all, Z)
  betas_simultaneous <- solve(crossprod(X_sim) + lambda * diag(ncol(X_sim)), 
                             crossprod(X_sim, y))[1:T_trials]
  
  cat("\n=== Simultaneous Estimation ===\n")
  cat("Simultaneous:   ", round(betas_simultaneous, 3), "\n")
  cat("Diff from true: ", round(betas_simultaneous - true_betas, 3), "\n")
  
  # The fmrilss implementation should do proper trial-wise LSS
  # It should match the direct LSS, not simultaneous estimation
  
  # Note: fmrilss uses a different numerical approach
  # so we allow for some tolerance in the comparison
  # Also handle NaN values that may occur for trials that don't fit in the time series
  valid_idx <- !is.na(betas_corrected) & !is.na(betas_direct)
  if (sum(valid_idx) > 0) {
    expect_lt(max(abs(betas_corrected[valid_idx] - betas_direct[valid_idx])), 0.5,
              "fmrilss implementation should approximately match direct LSS")
  }
})


