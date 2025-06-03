library(testthat)

# Helper implementing Least Squares Separate manually
manual_lss_separate <- function(y_proj, X_trials, h, P, lambda = 1e-6) {
  T_trials <- length(X_trials)
  betas <- numeric(T_trials)
  for (t in seq_len(T_trials)) {
    c_t <- P %*% (X_trials[[t]] %*% h)
    if (T_trials > 1) {
      c_other <- Reduce("+", lapply(X_trials[-t], function(X) P %*% (X %*% h)))
    } else {
      c_other <- matrix(0, nrow(P), 1)
    }
    X_full <- cbind(c_t, c_other)
    XtX <- crossprod(X_full) + lambda * diag(2)
    betas[t] <- solve(XtX, crossprod(X_full, y_proj))[1]
  }
  betas
}

test_that("run_lss_woodbury_corrected matches LS-S", {
  set.seed(42)
  n <- 60
  p <- 10
  T_trials <- 4

  X_trials <- lapply(seq_len(T_trials), function(t) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * (p + 2)
    if (onset + p - 1 <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })

  h <- exp(-(0:(p-1))/3)
  h <- h / sum(h)

  Z <- cbind(1, (1:n)/n)
  true_betas <- rnorm(T_trials)
  y <- Z %*% c(2, -1)
  for (t in seq_len(T_trials)) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.05)

  P <- diag(n) - Z %*% solve(crossprod(Z) + 1e-6 * diag(ncol(Z))) %*% t(Z)
  y_proj <- P %*% y

  beta_manual <- manual_lss_separate(y_proj, X_trials, h, P)
  beta_wood <- run_lss_woodbury_corrected(
    Y_proj_voxel_vector = as.vector(y_proj),
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    P_confound = P,
    lambda_ridge = 1e-6
  )
  expect_equal(beta_wood, beta_manual, tolerance = 1e-5)
})

test_that("run_lss_for_voxel_corrected_full matches LS-S", {
  set.seed(99)
  n <- 60
  p <- 10
  T_trials <- 3

  X_trials <- lapply(seq_len(T_trials), function(t) {
    X <- matrix(0, n, p)
    onset <- 6 + (t-1) * (p + 3)
    if (onset + p - 1 <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })

  h <- exp(-(0:(p-1))/3)
  h <- h / sum(h)

  Z <- cbind(1, (1:n)/n)
  true_betas <- rnorm(T_trials)
  y <- Z %*% c(1, -0.5)
  for (t in seq_len(T_trials)) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.05)

  P <- diag(n) - Z %*% solve(crossprod(Z) + 1e-6 * diag(ncol(Z))) %*% t(Z)
  y_proj <- P %*% y

  beta_manual <- manual_lss_separate(y_proj, X_trials, h, P)
  lss_prep <- prepare_lss_fixed_components_core(Z, 1, 1e-6)
  beta_full <- run_lss_for_voxel_corrected_full(
    Y_proj_voxel_vector = as.vector(y_proj),
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    A_lss_fixed_matrix = Z,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector
  )
  expect_equal(beta_full, beta_manual, tolerance = 1e-5)
})
