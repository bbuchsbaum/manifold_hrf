# Tests for the new unified interface

test_that("manifold_control works with different presets", {
  # Test balanced preset (default)
  ctrl_balanced <- manifold_control()
  expect_s3_class(ctrl_balanced, "manifold_control")
  expect_equal(attr(ctrl_balanced, "preset"), "balanced")
  expect_equal(ctrl_balanced$max_iter, 100)
  expect_equal(ctrl_balanced$num_neighbors_mfd, 75)
  
  # Test fast preset
  ctrl_fast <- manifold_control(preset = "fast")
  expect_equal(attr(ctrl_fast, "preset"), "fast")
  expect_equal(ctrl_fast$max_iter, 50)
  expect_equal(ctrl_fast$num_neighbors_mfd, 50)
  
  # Test thorough preset
  ctrl_thorough <- manifold_control(preset = "thorough")
  expect_equal(attr(ctrl_thorough, "preset"), "thorough")
  expect_equal(ctrl_thorough$max_iter, 500)
  expect_equal(ctrl_thorough$num_neighbors_mfd, 100)
})

test_that("manifold_control allows parameter overrides", {
  # Override specific parameters
  ctrl <- manifold_control(
    preset = "balanced",
    max_iter = 200,
    lambda_spatial_smooth = 1.0,
    verbose_level = 0
  )
  
  expect_equal(ctrl$max_iter, 200)  # Override
  expect_equal(ctrl$lambda_spatial_smooth, 1.0)  # Override
  expect_equal(ctrl$verbose_level, 0)  # Override
  expect_equal(ctrl$num_neighbors_mfd, 75)  # From balanced preset
})

test_that("manifold_control validates parameters", {
  # Invalid m_manifold_dim
  expect_error(
    manifold_control(m_manifold_dim = -1),
    "m_manifold_dim must be a positive integer"
  )
  
  # Invalid max_iter
  expect_error(
    manifold_control(max_iter = 0),
    "max_iter must be positive"
  )
  
  # Invalid lambda
  expect_error(
    manifold_control(lambda_manifold = -0.1),
    "lambda_manifold must be non-negative"
  )
  
  # Invalid neighbors
  expect_error(
    manifold_control(num_neighbors_mfd = 0),
    "num_neighbors_mfd must be at least 1"
  )
})

test_that("manifold_hrf_fit object creation works", {
  # Create test data
  set.seed(123)
  n <- 100
  V <- 10
  k <- 2
  m <- 3
  p <- 20
  
  # Create a test object
  fit <- new_manifold_hrf_fit(
    amplitudes = matrix(rnorm(k * V), k, V),
    trial_amplitudes = matrix(rnorm(50 * V), 50, V),
    hrf_shapes = matrix(rnorm(p * V), p, V),
    fitted_values = matrix(rnorm(n * V), n, V),
    residuals = matrix(rnorm(n * V), n, V),
    manifold_coords = matrix(rnorm(m * V), m, V),
    manifold = list(type = "test"),
    amplitudes_initial = matrix(rnorm(k * V), k, V),
    spatial_laplacian = diag(V),
    convergence_info = list(converged = TRUE, iterations = 10),
    call = quote(test_call()),
    control = manifold_control(),
    data_info = list(
      n_timepoints = n,
      n_voxels = V,
      n_trials = 50,
      n_conditions = k,
      TR = 2
    )
  )
  
  expect_s3_class(fit, "manifold_hrf_fit")
  expect_true(is.manifold_hrf_fit(fit))
})

test_that("manifold_hrf_fit methods work", {
  # Create a simple test object
  set.seed(456)
  n <- 50
  V <- 5
  k <- 2
  m <- 3
  p <- 10
  
  fit <- new_manifold_hrf_fit(
    amplitudes = matrix(1:10, k, V),
    trial_amplitudes = matrix(rnorm(20 * V), 20, V),
    hrf_shapes = matrix(rnorm(p * V), p, V),
    fitted_values = matrix(rnorm(n * V), n, V),
    residuals = matrix(rnorm(n * V, sd = 0.1), n, V),
    manifold_coords = matrix(rnorm(m * V), m, V),
    manifold = list(type = "test", dim = m),
    amplitudes_initial = matrix(rnorm(k * V), k, V),
    spatial_laplacian = diag(V),
    convergence_info = list(converged = TRUE, iterations = 5),
    call = quote(test_call()),
    control = manifold_control(preset = "fast"),
    data_info = list(
      n_timepoints = n,
      n_voxels = V,
      n_trials = 20,
      n_conditions = k,
      TR = 2
    )
  )
  
  # Test print method
  expect_output(print(fit), "Manifold HRF Fit")
  expect_output(print(fit), "Voxels: 5")
  expect_output(print(fit), "iterations: 5")
  
  # Test summary method
  summ <- summary(fit)
  expect_s3_class(summ, "summary.manifold_hrf_fit")
  expect_output(print(summ), "Summary of Manifold HRF Fit")
  
  # Test coef method
  coef_amp <- coef(fit, type = "amplitudes")
  expect_equal(dim(coef_amp), c(k, V))
  expect_equal(coef_amp, fit$amplitudes)
  
  coef_trial <- coef(fit, type = "trial_amplitudes")
  expect_equal(dim(coef_trial), c(20, V))
  
  coef_mfd <- coef(fit, type = "manifold_coords")
  expect_equal(dim(coef_mfd), c(m, V))
  
  # Test fitted and residuals methods
  fitted_vals <- fitted(fit)
  expect_equal(dim(fitted_vals), c(n, V))
  
  resid_vals <- residuals(fit)
  expect_equal(dim(resid_vals), c(n, V))
})

test_that("estimate_hrf_manifold validates inputs", {
  # Missing required columns in events
  expect_error(
    estimate_hrf_manifold(
      fmri_data = matrix(1, 100, 50),
      events = data.frame(time = 1:10),
      TR = 2
    ),
    "missing required columns"
  )
  
  # Invalid TR
  expect_error(
    estimate_hrf_manifold(
      fmri_data = matrix(1, 100, 50),
      events = data.frame(onset = 1:10, condition = rep("A", 10)),
      TR = -1
    ),
    "TR must be a positive number"
  )
  
  # Invalid control
  expect_error(
    estimate_hrf_manifold(
      fmri_data = matrix(1, 100, 50),
      events = data.frame(onset = 1:10, condition = rep("A", 10)),
      TR = 2,
      control = list(max_iter = 100)
    ),
    "control must be a manifold_control object"
  )
})

test_that("process_fmri_input handles different input types", {
  # Test with matrix input
  mat_data <- matrix(rnorm(100 * 50), 100, 50)
  result <- process_fmri_input(mat_data)
  
  expect_equal(dim(result$data_matrix), c(100, 50))
  expect_equal(dim(result$voxel_coords), c(50, 1))
  
  # Test with 4D array
  arr_data <- array(rnorm(10 * 10 * 5 * 100), c(10, 10, 5, 100))
  result_arr <- process_fmri_input(arr_data)
  
  expect_equal(nrow(result_arr$data_matrix), 100)
  expect_true(ncol(result_arr$data_matrix) <= 500)  # Some voxels may be zero
})

test_that("print.manifold_control works", {
  ctrl <- manifold_control(preset = "fast", verbose_level = 2)
  
  output <- capture.output(print(ctrl))
  expect_true(any(grepl("preset: fast", output)))
  expect_true(any(grepl("Manifold construction:", output)))
  expect_true(any(grepl("verbose_level: 2", output)))
})