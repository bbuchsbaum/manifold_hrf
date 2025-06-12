# Memory Heuristic Improvements Demo
# 
# This demo shows how the refactored memory heuristic calculation
# is more robust, testable, and well-engineered than the original

library(manifoldhrf)

# Example data dimensions (typical fMRI experiment)
n_timepoints <- 200   # 200 TRs
n_voxels <- 50000     # 50k voxels
n_trials <- 100       # 100 trials

cat("=== MEMORY HEURISTIC DEMO ===\n")
cat(sprintf("Dataset: %d timepoints × %d voxels × %d trials\n\n", 
           n_timepoints, n_voxels, n_trials))

# 1. ROBUST STRATEGY COMPARISON
cat("1. Strategy Comparison (with 32GB system):\n")
cat("─────────────────────────────────────────\n")

for (strategy in c("fast", "balanced", "memory")) {
  recommendation <- calculate_optimal_ram_heuristic(
    n_timepoints = n_timepoints,
    n_voxels = n_voxels, 
    n_trials = n_trials,
    target_strategy = strategy,
    system_memory_gb = 32.0,  # Simulate 32GB system
    verbose = FALSE
  )
  
  cat(sprintf("%-9s strategy: %.2f GB\n", 
             toupper(strategy), recommendation))
}

cat("\n2. Scalability Analysis:\n")
cat("─────────────────────────\n")

# Test different dataset sizes
test_cases <- data.frame(
  name = c("Small", "Medium", "Large", "Very Large"),
  voxels = c(10000, 50000, 100000, 200000),
  trials = c(50, 100, 200, 500)
)

for (i in seq_len(nrow(test_cases))) {
  case <- test_cases[i, ]
  
  rec <- calculate_optimal_ram_heuristic(
    n_timepoints = 200,
    n_voxels = case$voxels,
    n_trials = case$trials, 
    target_strategy = "balanced",
    system_memory_gb = 16.0,
    verbose = FALSE
  )
  
  cat(sprintf("%-10s (%5dk voxels, %3d trials): %.2f GB\n",
             case$name, case$voxels/1000, case$trials, rec))
}

cat("\n3. Memory Detection Demo:\n")
cat("─────────────────────────\n")

# This demonstrates the robust memory detection
detected_memory <- manifoldhrf:::.detect_system_memory_gb(verbose = TRUE)

if (!is.null(detected_memory)) {
  cat(sprintf("✓ Successfully detected %.1f GB system memory\n", detected_memory))
} else {
  cat("⚠ Memory detection failed, using fallback\n")
}

cat("\n4. Input Validation Demo:\n")
cat("─────────────────────────\n")

# Show proper error handling
tryCatch({
  calculate_optimal_ram_heuristic(100, 1000, 50, "invalid_strategy")
}, error = function(e) {
  cat("✓ Properly caught invalid strategy:", e$message, "\n")
})

# Show it works with valid inputs
tryCatch({
  result <- calculate_optimal_ram_heuristic(100, 1000, 50, "balanced", 
                                          system_memory_gb = 8, verbose = FALSE)
  cat("✓ Valid input produces result:", result, "GB\n")
}, error = function(e) {
  cat("✗ Unexpected error:", e$message, "\n")
})

cat("\n=== KEY IMPROVEMENTS ===\n")
cat("✓ Separation of concerns: detection, calculation, and configuration\n")
cat("✓ Robust error handling with graceful fallbacks\n") 
cat("✓ Cross-platform memory detection with multiple methods\n")
cat("✓ Input validation and sanity checks\n")
cat("✓ Testable design with dependency injection\n")
cat("✓ Clear configuration management\n")
cat("✓ Comprehensive logging and diagnostics\n")
cat("✓ No more magic numbers or hard-coded platform calls\n") 