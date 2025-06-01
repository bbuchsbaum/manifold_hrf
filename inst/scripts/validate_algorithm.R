# Simple validation script for M-HRF-LSS
# Tests the algorithm works on simulated data

library(manifoldhrf)

cat("=== M-HRF-LSS Algorithm Validation ===\n\n")

# Create simple simulated data
set.seed(123)
n <- 200  # timepoints
V <- 50   # voxels  
k <- 3    # conditions
p <- 24   # HRF length
TR <- 2

cat("Creating simulated fMRI data:\n")
cat(sprintf("  Timepoints: %d\n", n))
cat(sprintf("  Voxels: %d\n", V))
cat(sprintf("  Conditions: %d\n", k))
cat(sprintf("  TR: %d seconds\n\n", TR))

# Generate event timings
n_events_per_cond <- 10
event_list <- data.frame()

for (cond in 1:k) {
  onsets <- sort(runif(n_events_per_cond, min = 10, max = n * TR - 30))
  event_list <- rbind(event_list, data.frame(
    condition = paste0("Cond", cond),
    onset = onsets,
    duration = rep(1, n_events_per_cond)
  ))
}

# Create true HRF (double gamma)
time_hrf <- seq(0, by = TR, length.out = p)
true_hrf <- dgamma(time_hrf, shape = 6, scale = 0.9) - 
            0.35 * dgamma(time_hrf, shape = 16, scale = 0.9)
true_hrf <- true_hrf / sum(abs(true_hrf))

# Create design matrices
X_condition_list <- list()
for (cond in 1:k) {
  X_cond <- matrix(0, n, p)
  cond_events <- event_list[event_list$condition == paste0("Cond", cond), ]
  
  for (j in 1:nrow(cond_events)) {
    onset_idx <- round(cond_events$onset[j] / TR) + 1
    if (onset_idx <= n - p + 1) {
      X_cond[onset_idx:(onset_idx + p - 1), ] <- 
        X_cond[onset_idx:(onset_idx + p - 1), ] + diag(p)
    }
  }
  X_condition_list[[cond]] <- X_cond
}

# Generate data with known amplitudes
true_betas <- matrix(c(1, 1.5, 0.8), k, V)  # Different amplitudes per condition
true_betas <- true_betas + matrix(rnorm(k * V, sd = 0.1), k, V)  # Add variability

Y_signal <- matrix(0, n, V)
for (cond in 1:k) {
  # Convolve design with HRF
  X_conv <- matrix(0, n, V)
  design_pad <- rbind(X_condition_list[[cond]], matrix(0, p-1, p))
  
  for (t in 1:n) {
    for (v in 1:V) {
      X_conv[t, v] <- sum(design_pad[t:(t+p-1), ] %*% true_hrf)
    }
  }
  
  Y_signal <- Y_signal + X_conv * matrix(true_betas[cond, ], n, V, byrow = TRUE)
}

# Add noise
noise_sd <- 0.5
Y_data <- Y_signal + matrix(rnorm(n * V, sd = noise_sd), n, V)

cat("Data generated with:\n")
cat(sprintf("  SNR ≈ %.1f\n", sd(Y_signal) / noise_sd))
cat(sprintf("  Mean signal amplitude: %.2f\n\n", mean(abs(Y_signal))))

# Run M-HRF-LSS
cat("Running M-HRF-LSS algorithm...\n\n")

# Step 1: Create HRF library
cat("Step 1: Creating HRF library...\n")

# Create HRF library manually
time_points <- seq(0, by = TR, length.out = p)
N_hrfs <- 40

# Create grid of HRF parameters
shapes1 <- seq(4, 8, length.out = 10)
scales1 <- seq(0.8, 1.2, length.out = 4)

L_library <- matrix(0, p, N_hrfs)
idx <- 1

for (shape in shapes1) {
  for (scale in scales1) {
    if (idx <= N_hrfs) {
      # Double gamma HRF
      hrf <- dgamma(time_points, shape = shape, scale = scale) - 
             0.35 * dgamma(time_points, shape = shape * 2.5, scale = scale)
      # Normalize
      L_library[, idx] <- hrf / sum(abs(hrf))
      idx <- idx + 1
    }
  }
}

cat(sprintf("  Created library with %d HRF variants\n", ncol(L_library)))

# Step 2: Create manifold
cat("\nStep 2: Constructing HRF manifold...\n")

# Create affinity matrix
N <- ncol(L_library)
distances <- as.matrix(dist(t(L_library)))

# Self-tuning local scaling
k_nn <- 7
sigma_local <- numeric(N)
for (i in 1:N) {
  sorted_dists <- sort(distances[i, ])
  sigma_local[i] <- sorted_dists[k_nn + 1]
}

# Compute affinity
S_affinity <- matrix(0, N, N)
for (i in 1:(N-1)) {
  for (j in (i+1):N) {
    S_affinity[i, j] <- exp(-distances[i, j]^2 / (sigma_local[i] * sigma_local[j]))
    S_affinity[j, i] <- S_affinity[i, j]
  }
}

# Row normalize to get Markov matrix
S_markov <- S_affinity / rowSums(S_affinity)

# Get manifold basis
manifold <- get_manifold_basis_reconstructor_core(
  S_markov_matrix = S_markov,
  L_library_matrix = L_library,
  m_manifold_dim_target = 5,
  m_manifold_dim_min_variance = 0.95
)

cat(sprintf("  Manifold dimensions: %d\n", manifold$m_final_dim))
cat(sprintf("  Variance explained: %.1f%%\n", 
            sum(manifold$eigenvalues_S_vector[2:(manifold$m_final_dim + 1)]) * 100))

# Step 3: Voxel-wise HRF estimation
cat("\nStep 3: Estimating voxel-wise HRFs...\n")

# Transform designs to manifold basis
XB_list <- transform_designs_to_manifold_basis_core(
  X_condition_list_proj_matrices = X_condition_list,
  B_reconstructor_matrix = manifold$B_reconstructor_matrix
)

# Solve for gamma coefficients
Gamma <- solve_glm_for_gamma_core(
  Z_list_of_matrices = XB_list,
  Y_proj_matrix = Y_data,
  lambda_gamma = 0.01
)

# Extract manifold coordinates and amplitudes
svd_result <- extract_xi_beta_raw_svd_core(
  Gamma_coeffs_matrix = Gamma,
  m_manifold_dim = manifold$m_final_dim,
  k_conditions = k
)

# Apply identifiability constraints
ident_result <- apply_intrinsic_identifiability_core(
  Xi_raw_matrix = svd_result$Xi_raw_matrix,
  Beta_raw_matrix = svd_result$Beta_raw_matrix,
  Gamma_raw_matrix = Gamma
)

cat("  Voxel-wise estimation complete\n")

# Step 4: Spatial smoothing (minimal for this test)
cat("\nStep 4: Applying spatial smoothing...\n")
Xi_smooth <- apply_spatial_smoothing_core(
  Xi_ident_matrix = ident_result$Xi_ident_matrix,
  L_sp_sparse_matrix = diag(V),  # No spatial structure
  lambda_spatial_smooth = 0.01
)

# Step 5: Reconstruct HRF shapes
cat("\nStep 5: Reconstructing HRF shapes...\n")
hrf_shapes <- reconstruct_hrf_shapes_core(
  Xi_smooth_matrix = Xi_smooth,
  B_reconstructor_matrix = manifold$B_reconstructor_matrix
)

# Evaluate results
cat("\n=== RESULTS ===\n")

# 1. HRF recovery
hrf_correlations <- cor(true_hrf, hrf_shapes)
cat(sprintf("\nHRF Recovery:\n"))
cat(sprintf("  Median correlation with true HRF: %.3f\n", median(hrf_correlations)))
cat(sprintf("  Mean correlation: %.3f\n", mean(hrf_correlations)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(hrf_correlations), max(hrf_correlations)))

# 2. Amplitude recovery
estimated_betas <- ident_result$Beta_ident_matrix
beta_correlation <- cor(as.vector(true_betas), as.vector(estimated_betas))
cat(sprintf("\nAmplitude Recovery:\n"))
cat(sprintf("  Correlation with true amplitudes: %.3f\n", beta_correlation))
cat(sprintf("  Mean absolute error: %.3f\n", mean(abs(true_betas - estimated_betas))))

# 3. Reconstruction quality
Y_reconstructed <- matrix(0, n, V)
for (cond in 1:k) {
  X_conv <- matrix(0, n, V)
  design_pad <- rbind(X_condition_list[[cond]], matrix(0, p-1, p))
  
  for (t in 1:n) {
    for (v in 1:V) {
      X_conv[t, v] <- sum(design_pad[t:(t+p-1), ] %*% hrf_shapes[, v])
    }
  }
  
  Y_reconstructed <- Y_reconstructed + X_conv * 
    matrix(estimated_betas[cond, ], n, V, byrow = TRUE)
}

residuals <- Y_data - Y_reconstructed
r_squared <- 1 - sum(residuals^2) / sum((Y_data - mean(Y_data))^2)

cat(sprintf("\nReconstruction Quality:\n"))
cat(sprintf("  Overall R²: %.3f\n", r_squared))
cat(sprintf("  RMSE: %.3f\n", sqrt(mean(residuals^2))))

# 4. Algorithm confidence
confidence_score <- mean(c(
  median(hrf_correlations),  # HRF recovery
  beta_correlation,          # Amplitude recovery  
  r_squared                  # Reconstruction quality
))

cat(sprintf("\n=== OVERALL ALGORITHM CONFIDENCE: %.1f%% ===\n", confidence_score * 100))

if (confidence_score > 0.8) {
  cat("✓ Algorithm is working well!\n")
} else if (confidence_score > 0.6) {
  cat("⚠ Algorithm is working but may need parameter tuning\n")
} else {
  cat("✗ Algorithm performance is poor - investigation needed\n")
}

# Visualize sample results
if (interactive()) {
  par(mfrow = c(2, 2))
  
  # Plot 1: True vs estimated HRF (median voxel)
  median_voxel <- which.min(abs(hrf_correlations - median(hrf_correlations)))[1]
  plot(time_hrf, true_hrf, type = "l", lwd = 2, 
       main = "HRF Recovery (Median Voxel)",
       xlab = "Time (s)", ylab = "HRF", col = "black")
  lines(time_hrf, hrf_shapes[, median_voxel], col = "red", lwd = 2)
  legend("topright", c("True", "Estimated"), col = c("black", "red"), lwd = 2)
  
  # Plot 2: HRF correlation distribution
  hist(hrf_correlations, breaks = 20, main = "HRF Recovery Distribution",
       xlab = "Correlation with True HRF", col = "lightblue")
  abline(v = median(hrf_correlations), col = "red", lwd = 2)
  
  # Plot 3: True vs estimated amplitudes
  plot(true_betas, estimated_betas, 
       main = "Amplitude Recovery",
       xlab = "True Beta", ylab = "Estimated Beta",
       pch = 19, col = rep(1:k, each = V))
  abline(0, 1, col = "red", lwd = 2)
  
  # Plot 4: Sample time series
  plot(1:n, Y_data[, 1], type = "l", col = "gray",
       main = "Sample Voxel Fit",
       xlab = "Time", ylab = "Signal")
  lines(1:n, Y_reconstructed[, 1], col = "red", lwd = 2)
  legend("topright", c("Data", "Fitted"), col = c("gray", "red"), lwd = c(1, 2))
}

# Save results summary
results_summary <- list(
  hrf_recovery = list(
    correlations = hrf_correlations,
    median_cor = median(hrf_correlations)
  ),
  amplitude_recovery = list(
    correlation = beta_correlation,
    mae = mean(abs(true_betas - estimated_betas))
  ),
  reconstruction = list(
    r_squared = r_squared,
    rmse = sqrt(mean(residuals^2))
  ),
  confidence_score = confidence_score,
  parameters_used = list(
    lambda_gamma = 0.01,
    m_manifold = manifold$m_final_dim,
    n_hrfs = ncol(L_library)
  )
)

cat("\nValidation complete!\n")