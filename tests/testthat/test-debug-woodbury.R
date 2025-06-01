# Debug script for Woodbury LSS implementation
library(testthat)

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

# Woodbury
beta_wood <- run_lss_woodbury_corrected(
  Y_proj_voxel_vector = as.vector(y_proj),
  X_trial_onset_list_of_matrices = X_trials,
  H_shape_voxel_vector = hrf,
  lambda_ridge = 1e-6
)[trial]

# Naive
X_other <- do.call(cbind, lapply(X_trials[-trial], function(X) X %*% hrf))
X_full <- cbind(X_trials[[trial]] %*% hrf, X_other, Z)
XtX <- crossprod(X_full) + diag(1e-6, ncol(X_full))
Xty <- crossprod(X_full, y)
beta_naive <- solve(XtX, Xty)[1]

cat("Woodbury beta:", beta_wood, "\n")
cat("Naive beta:", beta_naive, "\n")
cat("True beta:", true_betas[trial], "\n")
cat("Difference:", abs(beta_wood - beta_naive), "\n")

# Debug: Check intermediate values
C_v <- do.call(cbind, lapply(X_trials, function(X) X %*% hrf))
U_v <- lss_prep$P_lss_matrix %*% C_v
V_reg <- C_v - Z %*% U_v

cat("\nDebug info:\n")
cat("P_lss dimensions:", dim(lss_prep$P_lss_matrix), "\n")
cat("C_v dimensions:", dim(C_v), "\n")
cat("U_v dimensions:", dim(U_v), "\n")
cat("V_reg dimensions:", dim(V_reg), "\n")
cat("p_lss vector length:", length(lss_prep$p_lss_vector), "\n")