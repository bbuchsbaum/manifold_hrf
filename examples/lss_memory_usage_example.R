# Example: Using LSS with Intelligent Memory Management

library(manifoldhrf)

# Typical fMRI dataset parameters
n_timepoints <- 300
n_voxels <- 50000  
n_trials <- 200
p_hrf <- 25

# 1. Check memory requirements
cat("\n=== Memory Requirement Analysis ===\n")

# Calculate memory for full precomputation
bytes_per_double <- 8
full_memory_gb <- (n_timepoints * n_trials * n_voxels * bytes_per_double) / (1024^3)

cat(sprintf("Dataset: %d timepoints × %d trials × %d voxels\n", 
            n_timepoints, n_trials, n_voxels))
cat(sprintf("Full precomputation would need: %.1f GB\n", full_memory_gb))

# 2. Calculate optimal ram_heuristic based on system
cat("\n=== System-Aware Configuration ===\n")

# Get recommended values for different strategies
fast_limit <- calculate_optimal_ram_heuristic(
  n_timepoints, n_voxels, n_trials, 
  target_strategy = "fast"
)

balanced_limit <- calculate_optimal_ram_heuristic(
  n_timepoints, n_voxels, n_trials,
  target_strategy = "balanced"  
)

memory_limit <- calculate_optimal_ram_heuristic(
  n_timepoints, n_voxels, n_trials,
  target_strategy = "memory"
)

# 3. Example: Running LSS with different memory strategies
cat("\n=== Running LSS with Different Strategies ===\n")

# Simulate some data for demonstration
set.seed(123)
Y_data <- matrix(rnorm(n_timepoints * 100), n_timepoints, 100)  # Just 100 voxels for demo
H_shapes <- matrix(rnorm(p_hrf * 100), p_hrf, 100)

# Create trial design matrices
X_trials <- lapply(1:20, function(t) {  # Just 20 trials for demo
  X <- matrix(0, n_timepoints, p_hrf)
  onset <- sample(1:(n_timepoints - p_hrf), 1)
  X[onset:(onset + p_hrf - 1), ] <- diag(p_hrf)
  X
})

# Prepare fixed components
A_fixed <- cbind(1, rnorm(n_timepoints))
lss_prep <- prepare_lss_fixed_components_core(A_fixed, lambda_fixed = 0.1)

# Strategy 1: Memory-conservative (forces streaming)
cat("\nStrategy 1: Memory-conservative\n")
system.time({
  beta_memory <- run_lss_voxel_loop_core_improved(
    Y_proj_matrix = Y_data,
    X_trial_onset_list_of_matrices = X_trials,
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = A_fixed,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    n_jobs = 1,
    ram_heuristic_GB_for_Rt = 0.001  # Force streaming
  )
})

# Strategy 2: Balanced (likely chunking)
cat("\nStrategy 2: Balanced\n")
system.time({
  beta_balanced <- run_lss_voxel_loop_core_improved(
    Y_proj_matrix = Y_data,
    X_trial_onset_list_of_matrices = X_trials,
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = A_fixed,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    n_jobs = 1,
    ram_heuristic_GB_for_Rt = balanced_limit
  )
})

# Strategy 3: Speed-optimized (full precomputation if possible)
cat("\nStrategy 3: Speed-optimized\n")
system.time({
  beta_fast <- run_lss_voxel_loop_core_improved(
    Y_proj_matrix = Y_data,
    X_trial_onset_list_of_matrices = X_trials,
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = A_fixed,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    n_jobs = 1,
    ram_heuristic_GB_for_Rt = 100  # Allow full precomputation
  )
})

# 4. Decision heuristic example
cat("\n=== Recommended Heuristic Logic ===\n")

decide_ram_heuristic <- function(n_timepoints, n_voxels, n_trials) {
  # Calculate dataset memory footprint
  data_gb <- (n_timepoints * n_trials * n_voxels * 8) / (1024^3)
  
  # Get available system memory (simplified)
  if (Sys.info()["sysname"] == "Darwin") {
    system_gb <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / (1024^3)
  } else {
    system_gb <- 16  # Conservative default
  }
  
  # Decision logic
  if (data_gb < 1) {
    # Small dataset: always precompute
    recommendation <- data_gb * 1.5
    strategy <- "full precomputation"
  } else if (data_gb < system_gb * 0.25) {
    # Fits comfortably: precompute
    recommendation <- data_gb * 1.2
    strategy <- "full precomputation"
  } else if (data_gb < system_gb * 0.5) {
    # Moderate: use chunking
    recommendation <- system_gb * 0.25
    strategy <- "chunked processing"
  } else {
    # Large: stream to avoid memory issues
    recommendation <- min(2.0, system_gb * 0.1)
    strategy <- "streaming"
  }
  
  cat(sprintf(
    "Dataset needs %.1f GB for full precomputation\n", data_gb
  ))
  cat(sprintf(
    "System has ~%.0f GB total RAM\n", system_gb
  ))
  cat(sprintf(
    "Recommendation: Set ram_heuristic_GB_for_Rt = %.1f\n", recommendation
  ))
  cat(sprintf(
    "This will use %s strategy\n", strategy
  ))
  
  return(recommendation)
}

# Example usage
recommended <- decide_ram_heuristic(n_timepoints, n_voxels, n_trials)

# 5. Using with the high-level interface
cat("\n=== Integration with mhrf_lss ===\n")

# The high-level function could be updated to use this automatically:
# mhrf_lss(..., 
#          ram_heuristic_GB_for_Rt = "auto")  # Would call decide_ram_heuristic

# Or explicitly:
# mhrf_lss(...,
#          ram_heuristic_GB_for_Rt = recommended)

# 6. Memory monitoring during execution
cat("\n=== Memory Usage Patterns ===\n")

# Function to monitor memory usage
monitor_memory <- function() {
  if (Sys.info()["sysname"] != "Windows") {
    # Unix-like systems
    pid <- Sys.getpid()
    mem_cmd <- sprintf("ps -p %d -o rss=", pid)
    mem_kb <- as.numeric(system(mem_cmd, intern = TRUE))
    mem_gb <- mem_kb / (1024^2)
  } else {
    # Windows
    mem_gb <- memory.size() / 1024
  }
  return(mem_gb)
}

# Example of monitoring during execution
cat("Memory usage examples:\n")
cat(sprintf("  Before LSS: %.2f GB\n", monitor_memory()))

# Run a small LSS computation
beta_small <- run_lss_voxel_loop_core(
  Y_proj_matrix = Y_data[, 1:10],
  X_trial_onset_list_of_matrices = X_trials[1:5],
  H_shapes_allvox_matrix = H_shapes[, 1:10],
  A_lss_fixed_matrix = A_fixed,
  P_lss_matrix = lss_prep$P_lss_matrix,
  p_lss_vector = lss_prep$p_lss_vector,
  n_jobs = 1
)

cat(sprintf("  After small LSS: %.2f GB\n", monitor_memory()))
gc()
cat(sprintf("  After garbage collection: %.2f GB\n", monitor_memory()))