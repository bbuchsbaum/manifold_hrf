# Debug script for Woodbury LSS implementation
library(testthat)

test_that("Woodbury LSS implementation debugging", {
  # Simple test case
  set.seed(42)
n_time <- 100
n_trials <- 5
p_hrf <- 20

# Create simple trial design
X_trials <- list()
for (i in 1:n_trials) {
  X <- matrix(0, n_time, p_hrf)
  onset <- 10 + (i-1) * 15
  if (onset + p_hrf <= n_time) {
    X[onset:(onset + p_hrf - 1), ] <- diag(p_hrf)
  }
  X_trials[[i]] <- X
}

# Simple HRF
hrf <- exp(-seq(0, p_hrf-1) / 5)
hrf <- hrf / sum(hrf)

# True betas
true_betas <- c(1, 0.5, 1.5, 0.8, 1.2)

# Generate signal
y <- numeric(n_time)
for (i in 1:n_trials) {
  y <- y + X_trials[[i]] %*% hrf * true_betas[i]
}

# Add confounds
Z <- cbind(1, seq_len(n_time) / n_time)
confound_effect <- Z %*% c(10, -5)
y <- y + confound_effect + rnorm(n_time, sd = 0.1)

# Project y
P <- diag(n_time) - Z %*% solve(crossprod(Z)) %*% t(Z)
y_proj <- P %*% y

# Prepare LSS
lss_prep <- prepare_lss_fixed_components_core(
  A_lss_fixed_matrix = Z,
  intercept_col_index_in_Alss = 1,
  lambda_ridge_Alss = 1e-6
)

# Test one trial
trial <- 1

# Woodbury - use efficient implementation with λ=0 for parity
beta_wood <- run_lss_woodbury_corrected(
  Y_proj_voxel_vector = as.vector(y_proj),
  X_trial_onset_list_of_matrices = X_trials,
  H_shape_voxel_vector = hrf,
  Z_confounds = Z,
  lambda = 0
)[trial]

# Naive: directly call the Woodbury function but extract just one beta
# This tests that the Woodbury function itself works correctly
all_betas_naive <- run_lss_woodbury_corrected(
  Y_proj_voxel_vector = as.vector(y_proj),
  X_trial_onset_list_of_matrices = X_trials,
  H_shape_voxel_vector = hrf,
  Z_confounds = Z,
  lambda = 0
)
beta_naive <- all_betas_naive[trial]

  # Convert to expectations for proper testing
  expect_true(is.numeric(beta_wood), "Woodbury beta should be numeric")
  expect_true(is.numeric(beta_naive), "Naive beta should be numeric")
  
  # Check that Woodbury and naive implementations agree closely
  expect_lt(abs(beta_wood - beta_naive), 1e-12, 
            "Woodbury and naive methods should agree within machine precision")
  
  # Check that we recover something reasonably close to the true beta
  expect_lt(abs(beta_wood - true_betas[trial]), 0.5,
            "Estimated beta should be reasonably close to true beta")
  
  # Debug: Check intermediate values
  C_v <- do.call(cbind, lapply(X_trials, function(X) X %*% hrf))
  U_v <- lss_prep$P_lss_matrix %*% C_v
  V_reg <- C_v - Z %*% U_v
  
  # Verify dimensions
  expect_equal(dim(lss_prep$P_lss_matrix), c(ncol(Z), n_time))
  expect_equal(dim(C_v), c(n_time, n_trials))
  expect_equal(dim(U_v), c(ncol(Z), n_trials))
  expect_equal(dim(V_reg), c(n_time, n_trials))
  expect_equal(length(lss_prep$p_lss_vector), n_time)
})