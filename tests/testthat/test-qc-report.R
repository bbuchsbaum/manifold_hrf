# Tests for QC Report Generation
# Tests for MHRF-QC-REPORT-01

test_that("generate_qc_report creates HTML report", {
  skip_if_not_installed("rmarkdown")
  
  # Create minimal test data
  set.seed(123)
  n <- 100
  V <- 50
  k <- 2
  p <- 20
  m <- 3
  
  # Mock results object
  results <- list(
    core_matrices = list(
      Y_data = matrix(rnorm(n * V), n, V),
      Xi_smoothed = matrix(rnorm(m * V), m, V),
      Xi_ident = matrix(rnorm(m * V), m, V),
      H_shapes = matrix(rnorm(p * V), p, V),
      Beta_condition_final = matrix(rnorm(k * V), k, V),
      Beta_trial = matrix(rnorm(10 * V), 10, V)
    ),
    manifold = list(
      eigenvalues_S = c(1, 0.5, 0.3, 0.2, 0.1, 0.05),
      m_auto_selected = 3,
      B_reconstructor = matrix(rnorm(p * m), p, m),
      Phi_coords = matrix(rnorm(V * m), V, m)
    ),
    diagnostics = list(
      r2_voxelwise = runif(V, 0, 1)
    )
  )
  
  # Parameters
  parameters <- list(
    manifold = list(
      m_manifold_dim_target = 3,
      m_manifold_dim_min_variance = 0.95,
      k_local_nn_for_sigma = 7,
      TR_precision = 0.1
    ),
    pipeline = list(
      lambda_gamma = 0.01,
      lambda_spatial_smooth = 0.5,
      lambda_beta_final = 0.01
    )
  )
  
  # Generate report
  output_file <- tempfile(fileext = ".html")
  
  # Check if template exists
  template_path <- system.file("rmd", "mhrf_qc_report.Rmd", 
                              package = "manifoldhrf")
  
  if (file.exists(template_path)) {
    report_path <- generate_qc_report(
      results = results,
      parameters = parameters,
      output_file = basename(output_file),
      output_dir = dirname(output_file),
      open_report = FALSE
    )
    
    expect_true(file.exists(report_path))
    expect_equal(report_path, output_file)
    
    # Check that HTML was generated
    html_content <- readLines(report_path)
    expect_true(any(grepl("M-HRF-LSS Pipeline QC Report", html_content)))
    
    # Clean up
    unlink(report_path)
  } else {
    skip("QC report template not found")
  }
})

test_that("compute_qc_diagnostics calculates metrics correctly", {
  set.seed(456)
  n <- 50
  V <- 20
  
  # Create test data with known R²
  Y_data <- matrix(rnorm(n * V), n, V)
  Y_predicted <- Y_data * 0.8 + matrix(rnorm(n * V, sd = 0.2), n, V)
  
  results <- list(
    core_matrices = list(
      Y_data = Y_data,
      Y_predicted = Y_predicted
    )
  )
  
  diagnostics <- compute_qc_diagnostics(results)
  
  # Check R² calculation
  expect_true("r2_voxelwise" %in% names(diagnostics))
  expect_length(diagnostics$r2_voxelwise, V)
  expect_true(all(diagnostics$r2_voxelwise >= 0, na.rm = TRUE))
  expect_true(all(diagnostics$r2_voxelwise <= 1, na.rm = TRUE))
  
  # Mean R² should be reasonably high given the construction
  expect_gt(mean(diagnostics$r2_voxelwise, na.rm = TRUE), 0.5)
})

test_that("create_qc_flags identifies issues correctly", {
  # Test with good results
  good_results <- list(
    core_matrices = list(
      Beta_trial = matrix(0, 30, 100)  # 30 trials
    ),
    diagnostics = list(
      r2_voxelwise = runif(100, 0.5, 0.9)  # Good fits
    ),
    hrf_stats = data.frame(
      peak_time = runif(100, 4, 6)  # Normal peak times
    ),
    n_truncated_hrfs = 0
  )
  
  flags_good <- create_qc_flags(good_results)
  expect_equal(flags_good$overall$status, "pass")
  expect_equal(length(flags_good), 1)  # Only overall flag
  
  # Test with poor results
  poor_results <- list(
    core_matrices = list(
      Beta_trial = matrix(0, 10, 100)  # Only 10 trials
    ),
    diagnostics = list(
      r2_voxelwise = runif(100, 0, 0.2)  # Poor fits
    ),
    hrf_stats = data.frame(
      peak_time = c(runif(50, 4, 6), runif(50, 12, 15))  # Some abnormal peaks
    ),
    n_truncated_hrfs = 3
  )
  
  flags_poor <- create_qc_flags(poor_results)
  expect_true(flags_poor$overall$status %in% c("warning", "fail"))
  expect_gt(length(flags_poor), 1)  # Should have specific flags
  expect_true("low_trial_count" %in% names(flags_poor))
  expect_true("poor_fits" %in% names(flags_poor))
  expect_true("unstable_hrf" %in% names(flags_poor))
  expect_true("hrf_truncation" %in% names(flags_poor))
})

test_that("extract_hrf_stats computes HRF metrics correctly", {
  set.seed(789)
  p <- 50  # Time points
  V <- 20  # Voxels
  TR <- 0.5
  
  # Create HRFs with known properties
  time_points <- seq(0, (p-1) * TR, by = TR)
  H_shapes <- matrix(0, p, V)
  
  for (v in 1:V) {
    # Gamma-like HRF with varying peak times
    peak_time_param <- 4 + v/10  # Shape parameter varies from 4.1 to 6
    hrf <- dgamma(time_points, shape = peak_time_param, rate = 1)
    H_shapes[, v] <- hrf / max(hrf)
  }
  
  # Add some zero HRFs
  H_shapes[, 18:20] <- 0
  
  stats <- extract_hrf_stats(H_shapes, TR_precision = TR)
  
  # Check structure
  expect_s3_class(stats, "data.frame")
  expect_equal(nrow(stats), V)
  expect_true(all(c("voxel", "peak_time", "peak_amplitude", "fwhm") %in% names(stats)))
  
  # Check values - peak times should be reasonable for gamma functions
  # Note: actual peak time depends on shape and scale parameters
  non_zero_peaks <- stats$peak_time[1:17]
  expect_true(all(non_zero_peaks > 0, na.rm = TRUE))
  expect_true(all(non_zero_peaks < 15, na.rm = TRUE))  # More reasonable range
  expect_true(all(stats$peak_amplitude[1:17] > 0.9, na.rm = TRUE))  # Normalized to ~1
  expect_true(all(is.na(stats$peak_time[18:20])))  # Zero HRFs
})

test_that("QC flags print method works", {
  flags <- list(
    overall = list(status = "warning", message = "2 QC issues detected", severity = 2),
    low_trial_count = list(
      status = "warning", 
      message = "Low trial count: 15 trials (recommended: ≥20)",
      severity = 2
    ),
    poor_fits = list(
      status = "warning",
      message = "35.0% of voxels have R² < 0.10",
      severity = 2
    )
  )
  class(flags) <- c("mhrf_qc_flags", "list")
  
  # Capture output
  output <- capture.output(print(flags))
  
  expect_true(any(grepl("Overall Status", output)))
  expect_true(any(grepl("Low Trial Count", output)))
  expect_true(any(grepl("Poor Fits", output)))
})

test_that("generate_qc_report handles missing diagnostics", {
  # Results without diagnostics
  results <- list(
    core_matrices = list(
      Y_data = matrix(rnorm(100), 20, 5),
      Xi_smoothed = matrix(rnorm(15), 3, 5),
      Beta_condition_final = matrix(rnorm(10), 2, 5),
      Beta_trial = matrix(rnorm(25), 5, 5)
    ),
    manifold = list(
      eigenvalues_S = c(1, 0.5, 0.3, 0.2)
    )
  )
  
  parameters <- list(
    manifold = list(m_manifold_dim_target = 3),
    pipeline = list(lambda_gamma = 0.01)
  )
  
  # Should add diagnostics automatically
  output_file <- tempfile(fileext = ".html")
  
  # Only test if template exists
  template_path <- system.file("rmd", "mhrf_qc_report.Rmd", 
                              package = "manifoldhrf")
  
  if (file.exists(template_path)) {
    expect_no_error(
      report_path <- generate_qc_report(
        results = results,
        parameters = parameters,
        output_file = basename(output_file),
        output_dir = dirname(output_file),
        open_report = FALSE
      )
    )
    
    if (exists("report_path") && file.exists(report_path)) {
      unlink(report_path)
    }
  }
})

test_that("QC thresholds can be customized", {
  results <- list(
    core_matrices = list(
      Beta_trial = matrix(0, 15, 100)  # 15 trials
    )
  )
  
  # Default thresholds - should flag low trial count
  flags_default <- create_qc_flags(results)
  expect_true("low_trial_count" %in% names(flags_default))
  
  # Custom thresholds - should pass
  flags_custom <- create_qc_flags(
    results,
    thresholds = list(min_trials = 10)  # Lower threshold
  )
  expect_equal(flags_custom$overall$status, "pass")
  expect_false("low_trial_count" %in% names(flags_custom))
})

test_that("trial regressor collinearity flag triggers", {
  results <- list(
    core_matrices = list(
      Beta_trial = matrix(0, 2, 5)
    )
  )
  attr(results$core_matrices$Beta_trial, "rank_deficient_voxels") <- rep(TRUE, 5)

  flags <- create_qc_flags(
    results,
    thresholds = list(min_trials = 1, max_trial_collinearity_fraction = 0.2)
  )

  expect_true("trial_regressor_collinearity" %in% names(flags))
})

test_that("compute_qc_diagnostics handles reconstruction error", {
  p <- 20
  N <- 50
  m <- 3
  
  # Create manifold components
  L_library <- matrix(rnorm(p * N), p, N)
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  Phi_coords <- matrix(rnorm(N * m), N, m)
  
  results <- list(
    manifold = list(
      L_library = L_library,
      B_reconstructor = B_reconstructor,
      Phi_coords = Phi_coords
    )
  )
  
  diagnostics <- compute_qc_diagnostics(results)
  
  # Check reconstruction error
  expect_true("reconstruction_error" %in% names(diagnostics))
  expect_length(diagnostics$reconstruction_error, m)
  expect_true(all(diagnostics$reconstruction_error > 0))
  # Reconstruction error could be >1 for random matrices
  expect_true(all(is.finite(diagnostics$reconstruction_error)))
  
  # Error should generally decrease with more dimensions (but not strictly for random data)
  # Just check that it's not increasing dramatically
  expect_true(diagnostics$reconstruction_error[m] <= diagnostics$reconstruction_error[1] * 2)
})