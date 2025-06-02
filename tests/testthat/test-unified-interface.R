# Integration tests for the unified mhrf_analyze interface

test_that("Basic mhrf_analyze interface works", {
  
  # Create minimal test data
  set.seed(42)
  n <- 100
  V <- 20
  Y_data <- matrix(rnorm(n * V), n, V)
  
  events <- data.frame(
    condition = rep(c("A", "B"), each = 5),
    onset = sort(runif(10, min = 10, max = 80)),
    duration = rep(1, 10)
  )
  
  # Test basic usage
  result <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = 2,
    verbose = 0
  )
  
  # Check result structure
  expect_s3_class(result, "mhrf_result")
  expect_true(is.list(result))
  
  # Check required components
  expect_true("hrf_shapes" %in% names(result))
  expect_true("amplitudes" %in% names(result))
  expect_true("manifold_coords" %in% names(result))
  expect_true("qc_metrics" %in% names(result))
  expect_true("metadata" %in% names(result))
  
  # Check dimensions
  expect_equal(ncol(result$hrf_shapes), V)
  expect_equal(ncol(result$amplitudes), V)
  expect_equal(nrow(result$amplitudes), 2)  # 2 conditions
  expect_equal(ncol(result$manifold_coords), V)
})


test_that("mhrf_analyze works with different presets", {
  
  set.seed(123)
  n <- 80
  V <- 15
  Y_data <- matrix(rnorm(n * V), n, V)
  
  events <- data.frame(
    condition = "task",
    onset = c(20, 40, 60),
    duration = 2
  )
  
  presets <- c("conservative", "balanced", "fast")
  
  for (preset in presets) {
    result <- mhrf_analyze(
      Y_data = Y_data,
      events = events,
      TR = 2,
      preset = preset,
      verbose = 0
    )
    
    expect_s3_class(result, "mhrf_result")
    expect_equal(result$metadata$preset_used, preset)
    expect_true(result$metadata$n_voxels_analyzed <= V)
  }
})


test_that("mhrf_analyze handles custom parameters", {
  
  set.seed(456)
  n <- 60
  V <- 10
  Y_data <- matrix(rnorm(n * V), n, V)
  
  events <- data.frame(
    condition = c("A", "B", "A"),
    onset = c(15, 30, 45),
    duration = 1
  )
  
  result <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = 2,
    lambda_gamma = 0.05,
    m_manifold_dim_target = 3,
    verbose = 0
  )
  
  expect_s3_class(result, "mhrf_result")
  expect_equal(result$metadata$parameters$lambda_gamma, 0.05)
  expect_equal(result$metadata$parameters$m_manifold_dim_target, 3)
})


test_that("mhrf_analyze S3 methods work correctly", {
  
  set.seed(789)
  n <- 50
  V <- 8
  Y_data <- matrix(rnorm(n * V), n, V)
  
  events <- data.frame(
    condition = c("rest", "task"),
    onset = c(10, 30),
    duration = c(5, 5)
  )
  
  result <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = 2,
    verbose = 0
  )
  
  # Test print method
  expect_output(print(result), "M-HRF-LSS Result")
  
  # Test summary method
  summary_result <- summary(result)
  expect_s3_class(summary_result, "summary.mhrf_result")
  expect_output(print(summary_result), "M-HRF-LSS Analysis Summary")
  
  # Test coef method
  amps <- coef(result, type = "amplitudes")
  expect_is(amps, "matrix")
  expect_equal(dim(amps), c(2, V))
  
  hrfs <- coef(result, type = "hrfs")
  expect_is(hrfs, "matrix")
  expect_equal(ncol(hrfs), V)
  
  # Test as.data.frame method
  df_amps <- as.data.frame(result, what = "amplitudes")
  expect_is(df_amps, "data.frame")
  expect_true("voxel" %in% names(df_amps))
  expect_true("condition" %in% names(df_amps))
  expect_true("amplitude" %in% names(df_amps))
  
  df_summary <- as.data.frame(result, what = "summary")
  expect_is(df_summary, "data.frame")
  expect_equal(nrow(df_summary), V)
})


test_that("mhrf_analyze error handling works", {
  
  # Test invalid Y_data
  expect_error(
    mhrf_analyze(
      Y_data = "not a matrix",
      events = data.frame(condition = "A", onset = 10, duration = 1),
      TR = 2
    ),
    "must be a matrix"
  )
  
  # Test missing required columns in events
  expect_error(
    mhrf_analyze(
      Y_data = matrix(rnorm(100), 50, 2),
      events = data.frame(time = 1:5),
      TR = 2
    ),
    "missing required columns"
  )
  
  # Test empty events
  expect_error(
    mhrf_analyze(
      Y_data = matrix(rnorm(100), 50, 2),
      events = data.frame(condition = character(0), onset = numeric(0), duration = numeric(0)),
      TR = 2
    )
  )
})


test_that("mhrf_analyze handles edge cases", {
  
  # Single voxel
  set.seed(111)
  Y_single <- matrix(rnorm(100), 100, 1)
  events_single <- data.frame(
    condition = "task",
    onset = c(20, 60),
    duration = 2
  )
  
  result_single <- mhrf_analyze(
    Y_data = Y_single,
    events = events_single,
    TR = 2,
    verbose = 0
  )
  
  expect_s3_class(result_single, "mhrf_result")
  expect_equal(ncol(result_single$hrf_shapes), 1)
  
  # Single condition
  Y_multi <- matrix(rnorm(200), 100, 2)
  events_single_cond <- data.frame(
    condition = "task",
    onset = c(10, 30, 50),
    duration = 1
  )
  
  result_single_cond <- mhrf_analyze(
    Y_data = Y_multi,
    events = events_single_cond,
    TR = 2,
    verbose = 0
  )
  
  expect_s3_class(result_single_cond, "mhrf_result")
  expect_equal(nrow(result_single_cond$amplitudes), 1)
})


test_that("mhrf_analyze trial-wise estimation works", {
  
  set.seed(222)
  n <- 120
  V <- 12
  Y_data <- matrix(rnorm(n * V), n, V)
  
  # Events with trial identifiers (sorted by onset)
  events_with_trials <- data.frame(
    condition = c("A", "B", "A", "B", "A", "B"),
    onset = c(10, 20, 30, 40, 50, 60),
    duration = 2,
    trial_type = paste0("trial_", 1:6)
  )
  
  result_trials <- mhrf_analyze(
    Y_data = Y_data,
    events = events_with_trials,
    TR = 2,
    verbose = 0
  )
  
  expect_s3_class(result_trials, "mhrf_result")
  expect_true(!is.null(result_trials$trial_amplitudes))
  expect_equal(result_trials$metadata$n_trials, 6)
})


test_that("mhrf_analyze respects voxel masking", {
  
  set.seed(333)
  n <- 80
  V <- 20
  Y_data <- matrix(rnorm(n * V), n, V)
  
  events <- data.frame(
    condition = "task",
    onset = c(15, 45),
    duration = 3
  )
  
  # Create mask that excludes half the voxels
  mask <- rep(c(TRUE, FALSE), V/2)
  
  result_masked <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = 2,
    voxel_mask = mask,
    verbose = 0
  )
  
  expect_s3_class(result_masked, "mhrf_result")
  expect_equal(result_masked$metadata$n_voxels_analyzed, sum(mask))
  expect_equal(result_masked$metadata$n_voxels_input, V)
})


test_that("mhrf_analyze progress tracking works", {
  
  set.seed(444)
  Y_data <- matrix(rnorm(400), 80, 5)
  events <- data.frame(
    condition = "A",
    onset = c(10, 40),
    duration = 2
  )
  
  # Test different verbosity levels
  expect_output(
    mhrf_analyze(Y_data, events, TR = 2, verbose = 1),
    "M-HRF-LSS Analysis"
  )
  
  # verbose = 0 should suppress progress messages but may still have warnings
  result <- expect_no_error({
    mhrf_analyze(Y_data, events, TR = 2, verbose = 0)
  })
  expect_type(result, "list")
})


test_that("mhrf_analyze metadata is complete", {
  
  set.seed(555)
  Y_data <- matrix(rnorm(300), 60, 5)
  events <- data.frame(
    condition = c("A", "B"),
    onset = c(15, 35),
    duration = 2
  )
  
  result <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = 2,
    preset = "balanced",
    verbose = 0
  )
  
  metadata <- result$metadata
  
  # Check all required metadata fields
  required_fields <- c(
    "n_timepoints", "n_voxels_input", "n_voxels_analyzed",
    "n_conditions", "n_trials", "n_hrfs_library",
    "manifold_dim", "manifold_method", "preset_used",
    "parameters", "runtime_seconds", "version"
  )
  
  for (field in required_fields) {
    expect_true(field %in% names(metadata), 
                info = paste("Missing metadata field:", field))
  }
  
  # Check metadata values
  expect_equal(metadata$n_timepoints, 60)
  expect_equal(metadata$n_voxels_input, 5)
  expect_equal(metadata$n_conditions, 2)
  expect_equal(metadata$preset_used, "balanced")
  expect_true(metadata$runtime_seconds >= 0)
})