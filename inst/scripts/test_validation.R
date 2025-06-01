# Test the improved input validation

library(manifoldhrf)

cat("=== Testing Input Validation ===\n\n")

# Test 1: Valid input (should work)
cat("Test 1: Valid input\n")
tryCatch({
  Y_valid <- matrix(rnorm(200), 100, 2)
  events_valid <- data.frame(
    condition = "task",
    onset = c(10, 40),
    duration = 2
  )
  
  result <- mhrf_analyze(Y_valid, events_valid, TR = 2, verbose = 0)
  cat("✓ Valid input accepted\n")
}, error = function(e) {
  cat("✗ Unexpected error with valid input:\n")
  cat(e$message, "\n")
})

# Test 2: Invalid Y_data type
cat("\nTest 2: Invalid Y_data type\n")
tryCatch({
  result <- mhrf_analyze("not a matrix", events_valid, TR = 2, verbose = 0)
  cat("✗ Should have failed\n")
}, error = function(e) {
  cat("✓ Correctly caught error:\n")
  cat(paste("  ", e$message), "\n")
})

# Test 3: Missing event columns
cat("\nTest 3: Missing event columns\n")
tryCatch({
  events_bad <- data.frame(time = c(10, 20), task = c("A", "B"))
  result <- mhrf_analyze(Y_valid, events_bad, TR = 2, verbose = 0)
  cat("✗ Should have failed\n")
}, error = function(e) {
  cat("✓ Correctly caught error:\n")
  cat(paste("  ", e$message), "\n")
})

# Test 4: Invalid TR
cat("\nTest 4: Invalid TR\n")
tryCatch({
  result <- mhrf_analyze(Y_valid, events_valid, TR = -1, verbose = 0)
  cat("✗ Should have failed\n")
}, error = function(e) {
  cat("✓ Correctly caught error:\n")
  cat(paste("  ", e$message), "\n")
})

# Test 5: Events after data end
cat("\nTest 5: Events after data end\n")
tryCatch({
  events_late <- data.frame(
    condition = "task",
    onset = c(10, 500),  # 500 seconds is way beyond 100*2=200 seconds of data
    duration = 2
  )
  result <- mhrf_analyze(Y_valid, events_late, TR = 2, verbose = 0)
  cat("✗ Should have failed\n")
}, error = function(e) {
  cat("✓ Correctly caught error:\n")
  cat(paste("  ", e$message), "\n")
})

# Test 6: Invalid preset
cat("\nTest 6: Invalid preset\n")
tryCatch({
  result <- mhrf_analyze(Y_valid, events_valid, TR = 2, preset = "invalid", verbose = 0)
  cat("✗ Should have failed\n")
}, error = function(e) {
  cat("✓ Correctly caught error:\n")
  cat(paste("  ", e$message), "\n")
})

# Test 7: Bad mask dimensions
cat("\nTest 7: Bad mask dimensions\n")
tryCatch({
  bad_mask <- c(TRUE, FALSE)  # Only 2 elements for 2 voxels data, should be OK
  result <- mhrf_analyze(Y_valid, events_valid, TR = 2, voxel_mask = bad_mask, verbose = 0)
  cat("✓ Good mask accepted\n")
}, error = function(e) {
  cat("Error with good mask:\n")
  cat(e$message, "\n")
})

tryCatch({
  bad_mask <- c(TRUE, FALSE, TRUE)  # 3 elements for 2 voxels - should fail
  result <- mhrf_analyze(Y_valid, events_valid, TR = 2, voxel_mask = bad_mask, verbose = 0)
  cat("✗ Should have failed\n")
}, error = function(e) {
  cat("✓ Correctly caught mask error:\n")
  cat(paste("  ", e$message), "\n")
})

# Test 8: Warning cases
cat("\nTest 8: Warning cases\n")

# Zero variance data
Y_zeros <- matrix(0, 50, 3)  # All zeros
tryCatch({
  result <- mhrf_analyze(Y_zeros, events_valid, TR = 2, verbose = 0)
  cat("✗ Should have failed with zero variance\n")
}, error = function(e) {
  cat("✓ Correctly caught zero variance error:\n")
  cat(paste("  ", e$message), "\n")
})

cat("\n=== Validation Testing Complete ===\n")