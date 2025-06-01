# Tests for Validation Simulation Framework
# Tests for MHRF-VALIDATE-SIM-01

test_that("run_mhrf_lss_simulation works with basic parameters", {
  # Run minimal simulation
  set.seed(42)
  
  sim_results <- run_mhrf_lss_simulation(
    n_voxels = 50,
    n_timepoints = 100,
    n_trials = 5,
    n_conditions = 2,
    TR = 2.0,
    noise_levels = c(0, 5),
    hrf_variability = "moderate",
    verbose = FALSE
  )
  
  # Check output structure
  expect_type(sim_results, "list")
  expect_true(all(c("ground_truth", "estimates", "metrics", 
                    "noise_curves", "report_path", "parameters") %in% names(sim_results)))
  
  # Check ground truth components
  expect_true(all(c("hrfs", "amplitudes", "design") %in% names(sim_results$ground_truth)))
  expect_equal(ncol(sim_results$ground_truth$hrfs$matrix), 50)  # n_voxels
  
  # Check metrics
  expect_s3_class(sim_results$metrics, "data.frame")
  expect_equal(nrow(sim_results$metrics), 2)  # One row per noise level
  expect_true("hrf_shape_correlation" %in% names(sim_results$metrics))
  expect_true("condition_amplitude_correlation" %in% names(sim_results$metrics))
  
  # Check parameters are stored correctly
  expect_equal(sim_results$parameters$simulation$n_voxels, 50)
  expect_equal(sim_results$parameters$simulation$noise_levels, c(0, 5))
})

test_that("generate_ground_truth_hrfs creates valid HRF library", {
  set.seed(123)
  
  # Test with no variability
  hrfs_none <- generate_ground_truth_hrfs(
    n_voxels = 20,
    hrf_variability = "none",
    TR = 2.0,
    manifold_params = list(TR_precision = 0.5)
  )
  
  expect_type(hrfs_none, "list")
  expect_equal(dim(hrfs_none$matrix)[2], 20)  # n_voxels
  expect_equal(length(hrfs_none$time_points), nrow(hrfs_none$matrix))
  expect_equal(hrfs_none$variability, "none")
  
  # Check HRFs are normalized
  max_vals <- apply(abs(hrfs_none$matrix), 2, max)
  expect_true(all(abs(max_vals - 1) < 0.01))
  
  # Test with high variability
  hrfs_high <- generate_ground_truth_hrfs(
    n_voxels = 30,
    hrf_variability = "high",
    TR = 1.0,
    manifold_params = list(TR_precision = 0.1)
  )
  
  # Check for more variation
  hrf_corrs <- cor(hrfs_high$matrix)
  diag(hrf_corrs) <- NA
  mean_corr <- mean(hrf_corrs, na.rm = TRUE)
  expect_lt(mean_corr, 0.98)  # Should have some variation (adjusted threshold)
})

test_that("generate_experimental_design creates valid designs", {
  design <- generate_experimental_design(
    n_timepoints = 200,
    n_trials = 10,
    n_conditions = 3,
    TR = 2.0
  )
  
  expect_type(design, "list")
  expect_equal(design$n_timepoints, 200)
  expect_equal(design$n_conditions, 3)
  expect_equal(length(design$conditions), design$total_trials)
  
  # Check design matrices
  expect_equal(length(design$X_condition_list), 3)
  expect_equal(length(design$X_trial_list), design$total_trials)
  
  # Check each condition matrix
  for (X in design$X_condition_list) {
    expect_equal(nrow(X), 200)
    expect_true(all(X >= 0))
    expect_true(sum(X) > 0)  # Should have some events
  }
  
  # Check onsets are reasonable
  expect_true(all(design$onsets > 0))
  expect_true(all(design$onsets < 200 * 2.0))  # Within experiment duration
})

test_that("generate_ground_truth_amplitudes creates realistic patterns", {
  amps <- generate_ground_truth_amplitudes(
    n_voxels = 60,
    n_conditions = 3,
    n_trials = 15,
    activation_patterns = c("sustained", "transient", "mixed")
  )
  
  expect_type(amps, "list")
  expect_equal(dim(amps$condition), c(3, 60))
  expect_equal(dim(amps$trial), c(15, 60))
  
  # Check that different patterns exist
  # Sustained pattern should have activity across conditions
  sustained_voxels <- 1:20
  expect_true(sum(abs(amps$condition[, sustained_voxels])) > 0)
  
  # Check trial amplitudes relate to condition amplitudes
  # (with some variability)
  cond1_active <- which(abs(amps$condition[1, ]) > 0.5)
  if (length(cond1_active) > 0) {
    # Trials for condition 1 should show activity in same voxels
    trial1 <- 1  # First trial (condition 1)
    expect_true(any(abs(amps$trial[trial1, cond1_active]) > 0))
  }
})

test_that("generate_bold_data creates realistic fMRI data", {
  set.seed(456)
  
  # Create simple test inputs
  hrfs <- list(
    matrix = matrix(1:20, nrow = 20, ncol = 10),
    time_points = seq(0, 19)
  )
  
  amps <- list(
    condition = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 10),
    trial = matrix(0.5, nrow = 6, ncol = 10)
  )
  
  design <- list(
    n_timepoints = 100,
    n_conditions = 2,
    X_condition_list = list(
      matrix(runif(100 * 20), 100, 20),
      matrix(runif(100 * 20), 100, 20)
    )
  )
  
  # Generate data with noise
  bold <- generate_bold_data(hrfs, amps, design, noise_level = 5, TR = 2.0)
  
  expect_type(bold, "list")
  expect_equal(dim(bold$Y_data), c(100, 10))
  expect_equal(dim(bold$Y_clean), c(100, 10))
  expect_true(ncol(bold$Z_confounds) >= 2)  # At least intercept and drift
  
  # Check DVARS is approximately correct
  dvars_actual <- compute_dvars(bold$Y_data)
  dvars_clean <- compute_dvars(bold$Y_clean)
  expect_gt(dvars_actual, dvars_clean)  # Noisy should have higher DVARS
})

test_that("compute_dvars calculates correctly", {
  # Simple test case
  Y <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  dvars <- compute_dvars(Y)
  
  expect_type(dvars, "double")
  expect_length(dvars, 1)
  expect_gt(dvars, 0)
  
  # Constant signal should have zero DVARS
  Y_const <- matrix(1, nrow = 10, ncol = 5)
  expect_equal(compute_dvars(Y_const), 0)
})

test_that("generate_fmri_noise creates temporally correlated noise", {
  noise <- generate_fmri_noise(
    n_timepoints = 100,
    n_voxels = 10,
    TR = 2.0
  )
  
  expect_equal(dim(noise), c(100, 10))
  
  # Check temporal autocorrelation
  # AR(1) noise should have positive lag-1 correlation
  lag1_cors <- sapply(1:10, function(v) {
    cor(noise[-100, v], noise[-1, v])
  })
  
  expect_true(mean(lag1_cors) > 0.1)  # Should have positive autocorrelation
  expect_true(mean(lag1_cors) < 0.5)  # But not too high
})

test_that("run_pipeline_on_simulated_data integrates all components", {
  skip_if_not_installed("Matrix")
  set.seed(789)
  
  # Create minimal test data
  n_voxels <- 20
  n_timepoints <- 80
  
  # Simple BOLD data
  bold_data <- list(
    Y_data = matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels),
    Z_confounds = cbind(1, seq_len(n_timepoints))
  )
  
  # Simple design
  design_info <- list(
    n_conditions = 2,
    X_condition_list = list(
      matrix(c(rep(1, 10), rep(0, 70)), n_timepoints, 10),
      matrix(c(rep(0, 40), rep(1, 10), rep(0, 30)), n_timepoints, 10)
    ),
    X_trial_list = list(
      matrix(c(rep(1, 10), rep(0, 70)), n_timepoints, 10),
      matrix(c(rep(0, 40), rep(1, 10), rep(0, 30)), n_timepoints, 10)
    )
  )
  
  # Simple HRFs
  ground_truth_hrfs <- list(
    matrix = matrix(rnorm(10 * n_voxels), 10, n_voxels)
  )
  
  # Run pipeline
  results <- run_pipeline_on_simulated_data(
    bold_data = bold_data,
    design_info = design_info,
    ground_truth_hrfs = ground_truth_hrfs,
    manifold_params = list(m_manifold_dim_target = 3),
    pipeline_params = list(lambda_gamma = 0.1),
    verbose = FALSE
  )
  
  expect_type(results, "list")
  expect_true("hrfs" %in% names(results))
  expect_true("beta_condition" %in% names(results))
  expect_true("beta_trial" %in% names(results))
  expect_equal(ncol(results$hrfs), n_voxels)
})

test_that("evaluate_pipeline_performance computes all metrics", {
  set.seed(321)
  
  # Create test estimates and ground truth
  n_voxels <- 30
  p <- 20
  k <- 2
  n_trials = 10
  
  estimates <- list(
    hrfs = matrix(rnorm(p * n_voxels), p, n_voxels),
    beta_condition = matrix(rnorm(k * n_voxels), k, n_voxels),
    beta_trial = matrix(rnorm(n_trials * n_voxels), n_trials, n_voxels),
    manifold_dim = 4
  )
  
  ground_truth_hrfs <- list(
    matrix = matrix(rnorm(p * n_voxels), p, n_voxels)
  )
  
  ground_truth_amplitudes <- list(
    condition = matrix(rnorm(k * n_voxels, mean = 1), k, n_voxels),
    trial = matrix(rnorm(n_trials * n_voxels, mean = 0.5), n_trials, n_voxels)
  )
  
  metrics <- evaluate_pipeline_performance(
    estimates, ground_truth_hrfs, ground_truth_amplitudes,
    noise_level = 2.5
  )
  
  expect_s3_class(metrics, "data.frame")
  expect_equal(nrow(metrics), 1)
  
  # Check key metrics exist
  expect_true("hrf_shape_correlation" %in% names(metrics))
  expect_true("condition_amplitude_correlation" %in% names(metrics))
  expect_true("trial_amplitude_correlation" %in% names(metrics))
  expect_true("spatial_dice_coefficient" %in% names(metrics))
  
  # Check metrics are in reasonable ranges
  expect_true(metrics$hrf_shape_correlation >= -1 && 
              metrics$hrf_shape_correlation <= 1)
  expect_true(metrics$spatial_dice_coefficient >= 0 && 
              metrics$spatial_dice_coefficient <= 1)
})

test_that("evaluate_hrf_recovery measures HRF similarity", {
  # Create test HRFs
  time_points <- seq(0, 20, by = 0.5)
  p <- length(time_points)
  n_voxels <- 10
  
  # True HRFs - gamma functions
  true_hrfs <- matrix(0, p, n_voxels)
  for (v in 1:n_voxels) {
    true_hrfs[, v] <- dgamma(time_points, shape = 6, rate = 1)
  }
  
  # Estimated HRFs - slightly shifted
  est_hrfs <- matrix(0, p, n_voxels)
  for (v in 1:n_voxels) {
    shift <- rnorm(1, 0, 0.5)
    est_hrfs[, v] <- dgamma(time_points - shift, shape = 6, rate = 1)
    est_hrfs[time_points < shift, v] <- 0
  }
  
  metrics <- evaluate_hrf_recovery(est_hrfs, true_hrfs)
  
  expect_type(metrics, "list")
  expect_true(all(c("peak_time_error", "fwhm_error", "shape_correlation") %in% 
                  names(metrics)))
  
  # Should have high correlation for similar shapes
  expect_gt(metrics$shape_correlation, 0.8)
  
  # Peak time error should be small
  expect_lt(metrics$peak_time_error, 2)
})

test_that("evaluate_amplitude_recovery measures beta accuracy", {
  # Perfect recovery case
  true_betas <- matrix(c(1, 0, 0, 1, 0.5, 0.5), nrow = 2, ncol = 3)
  est_betas <- true_betas + matrix(rnorm(6, sd = 0.1), 2, 3)
  
  metrics <- evaluate_amplitude_recovery(est_betas, true_betas, "condition")
  
  expect_type(metrics, "list")
  expect_gt(metrics$correlation, 0.9)
  expect_lt(metrics$rmse, 0.5)
  expect_gt(metrics$sensitivity, 0.5)
})

test_that("create_noise_robustness_curves organizes metrics", {
  # Create test metrics
  metrics_df <- data.frame(
    noise_level = c(0, 2, 5, 10),
    hrf_shape_correlation = c(0.95, 0.9, 0.8, 0.6),
    condition_amplitude_correlation = c(0.9, 0.85, 0.7, 0.5),
    trial_amplitude_correlation = c(0.85, 0.8, 0.65, 0.4),
    spatial_dice_coefficient = c(0.8, 0.75, 0.6, 0.4)
  )
  
  curves <- create_noise_robustness_curves(metrics_df)
  
  expect_type(curves, "list")
  expect_length(curves, 4)  # One for each metric
  
  # Check structure
  for (metric_name in names(curves)) {
    curve <- curves[[metric_name]]
    expect_true("data" %in% names(curve))
    expect_true("metric" %in% names(curve))
    expect_true("title" %in% names(curve))
    expect_equal(nrow(curve$data), 4)
  }
})

test_that("simulation handles edge cases gracefully", {
  # Very small simulation - may generate warnings about k neighbors
  small_sim <- suppressWarnings(
    run_mhrf_lss_simulation(
      n_voxels = 5,
      n_timepoints = 50,
      n_trials = 2,
      n_conditions = 2,
      noise_levels = c(0),
      verbose = FALSE
    )
  )
  
  expect_type(small_sim, "list")
  expect_true(file.exists(small_sim$report_path))
})

test_that("simulation parameters are validated", {
  # Invalid HRF variability
  expect_error(
    run_mhrf_lss_simulation(hrf_variability = "invalid"),
    "arg"
  )
  
  # Check other parameters work
  sim <- run_mhrf_lss_simulation(
    n_voxels = 20,
    n_timepoints = 60,
    noise_levels = c(1, 3),
    manifold_params = list(m_manifold_dim_target = 3),
    pipeline_params = list(lambda_gamma = 0.05),
    verbose = FALSE
  )
  
  expect_equal(sim$parameters$manifold$m_manifold_dim_target, 3)
  expect_equal(sim$parameters$pipeline$lambda_gamma, 0.05)
})

test_that("report generation works", {
  # Simple test data
  report_path <- generate_validation_report(
    ground_truth = list(hrfs = list(matrix = matrix(1, 10, 5))),
    all_results = list(noise_0 = list(hrfs = matrix(1, 10, 5))),
    metrics_df = data.frame(noise_level = 0, hrf_shape_correlation = 0.9),
    noise_curves = list(),
    params = list(),
    output_dir = tempdir()
  )
  
  expect_true(file.exists(report_path))
  expect_true(grepl("\\.rds$", report_path))
  
  # Load and check saved data
  saved_data <- readRDS(report_path)
  expect_type(saved_data, "list")
  expect_true("timestamp" %in% names(saved_data))
})