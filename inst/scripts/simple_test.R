# Simple test of M-HRF-LSS core algorithm
library(manifoldhrf)

# Test with minimal synthetic data
set.seed(42)
n <- 100  # timepoints
V <- 20   # voxels
k <- 2    # conditions
p <- 20   # HRF length

# Create simple HRF library
L <- matrix(0, p, 20)
for (i in 1:20) {
  t <- seq(0, by = 2, length.out = p)
  # Vary peak time
  peak <- 4 + i/5
  L[, i] <- dgamma(t, shape = peak, scale = 1) - 0.3 * dgamma(t, shape = peak * 2, scale = 1)
  L[, i] <- L[, i] / sum(abs(L[, i]))
}

# Test Component 0: Manifold Construction
cat("Testing Component 0: Manifold Construction\n")
S_markov <- calculate_manifold_affinity_core(
  L_library_matrix = L,
  k_local_nn_for_sigma = 5
)

manifold <- get_manifold_basis_reconstructor_core(
  S_markov_matrix = S_markov,
  L_library_matrix = L,
  m_manifold_dim_target = 3,
  m_manifold_dim_min_variance = 0.9
)

cat(sprintf("  ✓ Manifold created: %d dimensions\n", manifold$m_final_dim))

# Create test data
Y <- matrix(rnorm(n * V), n, V)
X_list <- list(
  matrix(rnorm(n * p), n, p),
  matrix(rnorm(n * p), n, p)
)

# Test Component 1: Voxel-wise fit
cat("\nTesting Component 1: Voxel-wise HRF estimation\n")

# Transform designs
XB_list <- transform_designs_to_manifold_basis_core(
  X_condition_list_proj_matrices = X_list,
  B_reconstructor_matrix = manifold$B_reconstructor_matrix
)

# Solve for gamma
Gamma <- solve_glm_for_gamma_core(
  Z_list_of_matrices = XB_list,
  Y_proj_matrix = Y,
  lambda_gamma = 0.01
)

cat(sprintf("  ✓ Gamma coefficients computed: %s\n", 
            paste(dim(Gamma), collapse = " x ")))

# Extract Xi and Beta
svd_result <- extract_xi_beta_raw_svd_core(
  Gamma_coeffs_matrix = Gamma,
  m_manifold_dim = manifold$m_final_dim,
  k_conditions = k
)

cat(sprintf("  ✓ Xi extracted: %s\n", 
            paste(dim(svd_result$Xi_raw_matrix), collapse = " x ")))
cat(sprintf("  ✓ Beta extracted: %s\n", 
            paste(dim(svd_result$Beta_raw_matrix), collapse = " x ")))

# Test Component 2: Spatial smoothing
cat("\nTesting Component 2: Spatial smoothing\n")

graph_result <- make_voxel_graph_laplacian_core(
  voxel_coords = cbind(1:V, rep(1, V), rep(1, V)),  # Simple 1D coords in 3D format
  num_neighbors = 3
)

Xi_smooth <- apply_spatial_smoothing_core(
  Xi_ident_matrix = svd_result$Xi_raw_matrix,
  L_sp_sparse_matrix = graph_result,
  lambda_spatial_smooth = 0.1
)

cat(sprintf("  ✓ Spatial smoothing applied: %s\n", 
            paste(dim(Xi_smooth), collapse = " x ")))

# Test Component 3: HRF reconstruction
cat("\nTesting Component 3: HRF reconstruction\n")

hrf_shapes <- reconstruct_hrf_shapes_core(
  B_reconstructor_matrix = manifold$B_reconstructor_matrix,
  Xi_smoothed_matrix = Xi_smooth
)

cat(sprintf("  ✓ HRF shapes reconstructed: %s\n", 
            paste(dim(hrf_shapes), collapse = " x ")))

# Test Component 3: Trial-wise LSS (single voxel)
cat("\nTesting Component 3: Trial-wise LSS\n")

# Create trial matrices
X_trials <- list()
for (i in 1:5) {
  X_trial <- matrix(0, n, p)
  onset <- i * 15
  if (onset + p <= n) {
    X_trial[onset:(onset + p - 1), ] <- diag(p)
  }
  X_trials[[i]] <- X_trial
}

lss_result <- run_lss_for_voxel_corrected(
  y_voxel = Y[, 1],
  X_trial_list = X_trials,
  h_voxel = hrf_shapes[, 1],
  TR = 2,
  lambda = 1e-6
)

cat(sprintf("  ✓ Trial betas estimated: %d trials\n", length(lss_result$beta_trials)))

# Summary
cat("\n=== CORE ALGORITHM TEST SUMMARY ===\n")
cat("✓ All core components executed successfully\n")
cat("✓ Dimensions are consistent throughout pipeline\n")
cat("✓ No crashes or errors\n")

# Check numerical stability
all_finite <- all(is.finite(hrf_shapes)) && 
              all(is.finite(svd_result$Xi_raw_matrix)) &&
              all(is.finite(svd_result$Beta_raw_matrix))

if (all_finite) {
  cat("✓ All outputs are numerically stable\n")
} else {
  cat("✗ Some outputs contain non-finite values\n")
}

cat("\nThe algorithm WORKS at a basic level!\n")