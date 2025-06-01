# Test the unified mhrf_analyze interface
library(manifoldhrf)

# Create simple simulated data
set.seed(42)
n <- 200  # timepoints
V <- 50   # voxels
TR <- 2

# Generate event timings
events <- data.frame(
  condition = rep(c("A", "B", "C"), each = 10),
  onset = sort(runif(30, min = 10, max = n * TR - 30)),
  duration = rep(1, 30)
)

# Generate simple data
Y_data <- matrix(rnorm(n * V), n, V)

# Add some signal
for (i in 1:nrow(events)) {
  onset_idx <- round(events$onset[i] / TR) + 1
  if (onset_idx <= n - 10) {
    cond_idx <- which(unique(events$condition) == events$condition[i])
    amplitude <- cond_idx * 0.5
    # Simple box-car response
    Y_data[onset_idx:(onset_idx + 5), ] <- Y_data[onset_idx:(onset_idx + 5), ] + amplitude
  }
}

cat("=== Testing Unified M-HRF-LSS Interface ===\n\n")

# Test 1: Basic usage
cat("Test 1: Basic usage with default parameters\n")
tryCatch({
  result1 <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = TR
  )
  cat("✓ Basic interface works\n")
  print(result1)
}, error = function(e) {
  cat("✗ Error in basic interface:\n")
  print(e)
})

# Test 2: Different presets
cat("\n\nTest 2: Testing different presets\n")
presets <- c("conservative", "balanced", "fast")

for (preset in presets) {
  cat(sprintf("\nTesting preset: %s\n", preset))
  tryCatch({
    result <- mhrf_analyze(
      Y_data = Y_data,
      events = events,
      TR = TR,
      preset = preset,
      verbose = 0  # Silent
    )
    cat(sprintf("✓ Preset '%s' works\n", preset))
  }, error = function(e) {
    cat(sprintf("✗ Error with preset '%s':\n", preset))
    print(e)
  })
}

# Test 3: Custom parameters
cat("\n\nTest 3: Testing custom parameters\n")
tryCatch({
  result3 <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = TR,
    preset = "balanced",
    lambda_gamma = 0.05,
    m_manifold_dim_target = 4,
    verbose = 1
  )
  cat("✓ Custom parameters work\n")
}, error = function(e) {
  cat("✗ Error with custom parameters:\n")
  print(e)
})

# Test 4: Result methods
cat("\n\nTest 4: Testing result S3 methods\n")
if (exists("result1")) {
  tryCatch({
    # Summary
    cat("\nSummary method:\n")
    summary(result1)
    
    # Coefficients
    cat("\nExtract coefficients:\n")
    amps <- coef(result1, type = "amplitudes")
    cat(sprintf("  Amplitudes dimensions: %s\n", paste(dim(amps), collapse = " x ")))
    
    # As data frame
    cat("\nConvert to data frame:\n")
    df <- as.data.frame(result1, what = "summary")
    cat(sprintf("  Data frame rows: %d\n", nrow(df)))
    
    cat("✓ All S3 methods work\n")
  }, error = function(e) {
    cat("✗ Error with S3 methods:\n")
    print(e)
  })
}

# Test 5: Error handling
cat("\n\nTest 5: Testing error handling\n")

# Bad input
tryCatch({
  result_bad <- mhrf_analyze(
    Y_data = "not a matrix",
    events = events,
    TR = TR
  )
}, error = function(e) {
  cat("✓ Correctly caught bad input error\n")
})

# Missing columns in events
tryCatch({
  bad_events <- data.frame(time = 1:10)
  result_bad <- mhrf_analyze(
    Y_data = Y_data,
    events = bad_events,
    TR = TR
  )
}, error = function(e) {
  cat("✓ Correctly caught bad events error\n")
})

cat("\n=== Interface Testing Complete ===\n")