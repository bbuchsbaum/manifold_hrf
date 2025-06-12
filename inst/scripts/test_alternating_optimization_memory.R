# Memory usage comparison for alternating optimization refactoring
library(manifoldhrf)

# Set up test parameters mimicking whole brain analysis
n <- 200   # timepoints (typical fMRI)
p <- 30    # HRF length
k <- 4     # conditions
V <- 50000 # voxels (whole brain)

cat("Test parameters:\n")
cat(sprintf("- Timepoints (n): %d\n", n))
cat(sprintf("- HRF length (p): %d\n", p))
cat(sprintf("- Conditions (k): %d\n", k))
cat(sprintf("- Voxels (V): %d\n", V))
cat("\n")

# Calculate memory requirements
bytes_per_double <- 8

# Old implementation would require:
old_memory_GB <- k * n * V * bytes_per_double / 1e9
cat(sprintf("Old implementation memory requirement: %.2f GB\n", old_memory_GB))
cat("  (Stores k × (n × V) matrices for all condition regressors)\n\n")

# New streaming implementation requires:
new_memory_GB <- n * k * bytes_per_double / 1e9  # Just one n×k matrix per voxel
cat(sprintf("New streaming implementation memory requirement: %.3f MB\n", new_memory_GB * 1000))
cat("  (Computes n × k design matrix on-the-fly per voxel)\n\n")

# Memory reduction factor
reduction_factor <- old_memory_GB / new_memory_GB
cat(sprintf("Memory reduction factor: %.0fx\n", reduction_factor))
cat(sprintf("Memory saved: %.2f GB\n", old_memory_GB - new_memory_GB))
cat("\n")

# Demonstrate with smaller example that actually runs
if (interactive()) {
  cat("Running small demonstration...\n")
  
  # Small demo parameters
  n_demo <- 100
  V_demo <- 500
  k_demo <- 3
  p_demo <- 20
  
  # Generate test data
  set.seed(123)
  Y_proj <- matrix(rnorm(n_demo * V_demo), n_demo, V_demo)
  X_cond_list <- lapply(1:k_demo, function(c) {
    matrix(rnorm(n_demo * p_demo) / 10, n_demo, p_demo)
  })
  H_shapes <- matrix(rnorm(p_demo * V_demo), p_demo, V_demo)
  
  # Time the new implementation
  gc(full = TRUE)
  time_start <- Sys.time()
  
  Beta_result <- estimate_final_condition_betas_core(
    Y_proj, X_cond_list, H_shapes,
    lambda_beta_final = 0.01,
    n_jobs = 1
  )
  
  time_end <- Sys.time()
  elapsed <- as.numeric(time_end - time_start, units = "secs")
  
  cat(sprintf("\nDemo completed in %.2f seconds\n", elapsed))
  cat(sprintf("Result dimensions: %d × %d\n", nrow(Beta_result), ncol(Beta_result)))
  cat(sprintf("All values finite: %s\n", all(is.finite(Beta_result))))
}