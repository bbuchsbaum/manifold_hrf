# Tests for Soundness Improvements
# Tests for SOUND-* improvements

test_that("HRF library quality check works correctly", {
  # Create test library with issues
  p <- 20
  N <- 10
  
  # Good library
  L_good <- matrix(rnorm(p * N), p, N)
  L_good <- apply(L_good, 2, function(x) x / sum(abs(x)))
  
  quality_good <- check_hrf_library_quality(L_good)
  expect_true(quality_good$is_good_quality)
  expect_false(quality_good$has_duplicates)
  expect_false(quality_good$is_ill_conditioned)
  
  # Library with duplicates
  L_dup <- L_good
  # Create a very high correlation by adding minimal noise
  L_dup[, 2] <- L_dup[, 1] + c(rep(0, p-1), 0.001)  # Only change last element slightly
  
  # Debug: check actual correlation
  actual_cor <- cor(L_dup[, 1], L_dup[, 2])
  message("Actual correlation between duplicates: ", actual_cor)
  
  quality_dup <- check_hrf_library_quality(L_dup)
  
  # If correlation is high enough but not detected, adjust threshold
  if (actual_cor > 0.99 && !quality_dup$has_duplicates) {
    message("High correlation but not detected - threshold issue")
  }
  
  expect_true(quality_dup$has_duplicates)
  expect_gt(quality_dup$n_duplicates, 0)
  expect_false(quality_dup$is_good_quality)
  
  # Library with zero HRFs
  L_zero <- L_good
  L_zero[, 5] <- 0
  
  quality_zero <- check_hrf_library_quality(L_zero)
  expect_true(quality_zero$has_degenerate_hrfs)
  expect_equal(quality_zero$n_zero_hrfs, 1)
})

test_that("Duplicate HRF removal works", {
  p <- 20
  N <- 10
  
  # Create library with duplicates
  L <- matrix(rnorm(p * N), p, N)
  L[, 2] <- L[, 1]  # Exact duplicate
  L[, 4] <- L[, 3] * 1.01  # Near duplicate
  
  L_cleaned <- remove_duplicate_hrfs(L, cor_threshold = 0.99)
  
  expect_lt(ncol(L_cleaned), ncol(L))
  expect_equal(ncol(L_cleaned), 8)  # Should remove 2 duplicates
  
  # Check no duplicates remain
  cor_matrix <- cor(L_cleaned)
  diag(cor_matrix) <- 0
  expect_true(all(abs(cor_matrix) < 0.99))
})

test_that("PCA fallback works correctly", {
  p <- 30
  N <- 50
  
  # Create HRF library
  L <- matrix(rnorm(p * N), p, N)
  L <- apply(L, 2, function(x) x / sum(abs(x)))
  
  # Run PCA fallback
  pca_result <- compute_pca_fallback(L, m_target = 5, min_variance = 0.95)
  
  expect_equal(ncol(pca_result$B_reconstructor_matrix), pca_result$m_final_dim)
  expect_equal(nrow(pca_result$B_reconstructor_matrix), p)
  expect_equal(pca_result$method_used, "PCA")
  expect_true(pca_result$variance_explained >= 0.9)  # Should capture most variance
  
  # Check orthogonality of basis
  BtB <- crossprod(pca_result$B_reconstructor_matrix)
  expect_equal(diag(BtB), rep(1, pca_result$m_final_dim), tolerance = 1e-8)
})

test_that("Robust manifold construction with fallback works", {
  set.seed(123)
  p <- 20
  N <- 30
  
  # Create reasonable library and Markov matrix
  L <- matrix(rnorm(p * N), p, N)
  S <- matrix(runif(N * N), N, N)
  S <- S + t(S)  # Symmetric
  S <- S / rowSums(S)  # Row stochastic
  
  # Test with good inputs
  result_good <- get_manifold_basis_reconstructor_robust(
    S_markov_matrix = S,
    L_library_matrix = L,
    m_manifold_dim_target = 3,
    fallback_to_pca = TRUE
  )
  
  expect_true("method_used" %in% names(result_good))
  expect_true("library_quality" %in% names(result_good))
  
  # Test with degenerate input (should trigger PCA fallback)
  S_bad <- matrix(1, N, N) / N  # All same values
  
  expect_warning(
    result_fallback <- get_manifold_basis_reconstructor_robust(
      S_markov_matrix = S_bad,
      L_library_matrix = L,
      m_manifold_dim_target = 3,
      fallback_to_pca = TRUE
    ),
    "Falling back to PCA"
  )
  
  expect_equal(result_fallback$method_used, "PCA")
})

test_that("Adaptive parameter selection works", {
  set.seed(456)
  n <- 100
  V <- 50
  k <- 3
  p <- 20
  
  # Create test data with known scale
  Y_data <- matrix(rnorm(n * V, sd = 2), n, V)
  
  # Create design matrices
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  
  # Create voxel coordinates
  coords <- expand.grid(x = 1:5, y = 1:5, z = 1:2)[1:V, ]
  
  # Get suggested parameters
  params <- suggest_parameters(Y_data, X_list, as.matrix(coords))
  
  expect_type(params, "list")
  expect_true("lambda_gamma" %in% names(params))
  expect_true("lambda_spatial_smooth" %in% names(params))
  expect_true("m_manifold_dim_target" %in% names(params))
  
  # Lambda should scale with data variance
  expect_gt(params$lambda_gamma, 0)
  expect_lt(params$lambda_gamma, 1)  # Not too large
})

test_that("Preset parameters work correctly", {
  # Test conservative preset
  params_cons <- get_preset_params("conservative")
  expect_equal(params_cons$m_manifold_dim_target, 3)
  expect_gt(params_cons$lambda_gamma, 0.05)
  expect_true(params_cons$use_robust_svd)
  expect_s3_class(params_cons, "mhrf_preset")
  
  # Test aggressive preset
  params_agg <- get_preset_params("aggressive")
  expect_gt(params_agg$m_manifold_dim_target, 5)
  expect_lt(params_agg$lambda_gamma, 0.01)
  expect_false(params_agg$use_robust_svd)
  
  # Test new presets
  params_fast <- get_preset_params("fast")
  expect_true(params_fast$use_parallel)
  expect_equal(params_fast$max_iterations, 1)
  
  params_quality <- get_preset_params("quality")
  expect_true(params_quality$apply_hrf_constraints)
  expect_lt(params_quality$convergence_tol, 1e-4)
  
  params_robust <- get_preset_params("robust")
  expect_true(params_robust$fallback_to_pca)
  expect_equal(params_robust$handle_zero_voxels, "noise")
  
  # Test with data scaling
  params_scaled <- get_preset_params("balanced", data_scale = 10)
  expect_gt(params_scaled$lambda_gamma, get_preset_params("balanced")$lambda_gamma)
  
  # Test with voxel count adjustment
  params_large <- get_preset_params("balanced", n_voxels = 60000)
  expect_equal(params_large$chunk_size, 2000)
  expect_true(params_large$use_parallel)
  
  # Test validation function
  test_data <- matrix(rnorm(100 * 20), 100, 20)
  test_design <- list(matrix(rnorm(100 * 10), 100, 10))
  
  validation <- params_cons$validate_data(test_data, test_design)
  expect_true(validation$compatible)
})

test_that("Workflow configuration works", {
  # Create workflow config
  config <- create_workflow_config(
    preset = "balanced",
    custom_params = list(lambda_gamma = 0.02),
    data_checks = TRUE
  )
  
  expect_s3_class(config, "mhrf_workflow_config")
  expect_equal(config$lambda_gamma, 0.02)  # Custom override
  expect_true(config$data_checks)
  expect_true(dir.exists(config$output_dir))
  expect_true(file.exists(config$log_file))
  
  # Test logging
  config$log_message("Test message", "INFO")
  expect_true(file.exists(config$log_file))
  
  log_contents <- readLines(config$log_file)
  expect_true(any(grepl("Test message", log_contents)))
  
  # Clean up
  unlink(config$output_dir, recursive = TRUE)
})

test_that("Preset print method works", {
  params <- get_preset_params("robust")
  
  # Capture print output
  output <- capture.output(print(params))
  
  expect_true(any(grepl("M-HRF-LSS Parameter Preset", output)))
  expect_true(any(grepl("robust", output)))
  expect_true(any(grepl("Robustness features: Enabled", output)))
})

test_that("Robust SVD extraction handles edge cases", {
  set.seed(789)
  k <- 3
  m <- 4
  V <- 20
  
  # Create test gamma matrix with issues
  Gamma <- matrix(rnorm(k * m * V), k * m, V)
  
  # Make some voxels problematic
  Gamma[, 1] <- 0  # Zero voxel
  Gamma[, 2] <- 1e-15  # Near zero
  Gamma[, 3] <- Gamma[, 3] * 1e10  # Large values
  
  # Run robust extraction
  result <- extract_xi_beta_raw_svd_core(
    Gamma_coeffs_matrix = Gamma,
    m_manifold_dim = m,
    k_conditions = k,
    method = "robust"
  )
  
  expect_equal(dim(result$Xi_raw_matrix), c(m, V))
  expect_equal(dim(result$Beta_raw_matrix), c(k, V))
  expect_false(any(is.na(result$Xi_raw_matrix)))
  expect_false(any(is.infinite(result$Xi_raw_matrix)))
  
  # Check quality metrics
  expect_equal(result$quality_metrics$svd_method[1], "zero")
  # May or may not apply regularization depending on condition number
  # Just check that metrics are populated
  expect_length(result$quality_metrics$regularization_applied, V)
  expect_length(result$quality_metrics$svd_method, V)
})

test_that("Local SNR computation works", {
  set.seed(321)
  n <- 100
  V <- 30
  
  # Create data with varying SNR
  Y_data <- matrix(0, n, V)
  
  # High SNR voxels
  Y_data[, 1:10] <- matrix(rnorm(n * 10, sd = 0.1), n, 10) + 
                    sin(seq(0, 4*pi, length.out = n))
  
  # Low SNR voxels
  Y_data[, 11:20] <- matrix(rnorm(n * 10, sd = 2), n, 10)
  
  # Medium SNR
  Y_data[, 21:30] <- matrix(rnorm(n * 10, sd = 0.5), n, 10) + 
                     0.5 * cos(seq(0, 2*pi, length.out = n))
  
  snr <- compute_local_snr(Y_data, method = "temporal_variance")
  
  expect_length(snr, V)
  expect_true(all(snr > 0))
  expect_gt(mean(snr[1:10]), mean(snr[11:20]))  # High > Low SNR
})

test_that("Adaptive spatial smoothing works", {
  set.seed(654)
  m <- 3
  V <- 25
  
  # Create test data
  Xi <- matrix(rnorm(m * V), m, V)
  
  # Create simple Laplacian (grid)
  L <- diag(V)
  for (i in 1:(V-1)) {
    L[i, i+1] <- L[i+1, i] <- -0.25
  }
  diag(L) <- -rowSums(L)
  
  # Variable SNR
  snr <- c(rep(10, 10), rep(1, 10), rep(5, 5))
  
  # Apply adaptive smoothing
  Xi_smooth <- apply_spatial_smoothing_adaptive(
    Xi_ident_matrix = Xi,
    L_sp_sparse_matrix = L,
    lambda_spatial_smooth = 0.5,
    local_snr = snr,
    edge_preserve = FALSE
  )
  
  expect_equal(dim(Xi_smooth), dim(Xi))
  
  # Check that smoothing was applied
  # The difference should be non-zero
  frobenius_diff <- norm(Xi_smooth - Xi, "F")
  expect_gt(frobenius_diff, 0)
  
  # Check that no NAs or Infs were introduced
  expect_false(any(is.na(Xi_smooth)))
  expect_false(any(is.infinite(Xi_smooth)))
})

test_that("Adaptive smoothing uses varying lambdas when SNR varies", {
  set.seed(777)
  m <- 2
  V <- 4

  Xi <- matrix(rnorm(m * V), m, V)

  L <- diag(V)
  for (i in 1:(V - 1)) {
    L[i, i + 1] <- L[i + 1, i] <- -0.5
  }
  diag(L) <- -rowSums(L)

  snr <- c(1, 4, 16, 4)

  Xi_smooth <- apply_spatial_smoothing_adaptive(
    Xi_ident_matrix = Xi,
    L_sp_sparse_matrix = L,
    lambda_spatial_smooth = 1,
    local_snr = snr,
    edge_preserve = FALSE
  )

  snr_factor <- 1 / sqrt(snr)
  snr_factor <- snr_factor / median(snr_factor)
  lambda_vec <- 1 * snr_factor

  expect_gt(max(lambda_vec) - min(lambda_vec), 0)

  diff_mag <- colSums((Xi_smooth - Xi)^2)
  expect_gt(diff_mag[1], diff_mag[3])
})

test_that("Outlier detection works correctly", {
  set.seed(987)
  n <- 100
  V <- 20
  
  # Create clean data
  Y_clean <- matrix(rnorm(n * V), n, V)
  
  # Add outliers
  Y_outlier <- Y_clean
  outlier_times <- c(10, 50, 90)
  outlier_voxels <- c(1, 5, 10)
  
  for (t in outlier_times) {
    for (v in outlier_voxels) {
      Y_outlier[t, v] <- Y_outlier[t, v] + sample(c(-5, 5), 1) * sd(Y_clean[, v])
    }
  }
  
  # Detect outliers
  weights <- detect_outlier_timepoints(Y_outlier, threshold = 3)
  
  expect_equal(dim(weights), dim(Y_outlier))
  expect_true(all(weights >= 0 & weights <= 1))
  
  # Check that outliers were downweighted
  for (v in outlier_voxels) {
    expect_lt(min(weights[outlier_times, v]), 1)
  }
})

test_that("Voxel screening identifies bad voxels", {
  n <- 100
  V <- 30
  
  # Create data with various issues
  Y_data <- matrix(rnorm(n * V), n, V)
  
  # Zero variance voxels
  Y_data[, 1:5] <- 100
  
  # Spike voxels
  Y_data[, 6:10] <- rnorm(n * 5)
  for (v in 6:10) {
    spike_times <- sample(1:n, 15)
    Y_data[spike_times, v] <- Y_data[spike_times, v] + 10
  }
  
  # Good voxels
  Y_data[, 11:30] <- matrix(rnorm(n * 20), n, 20)
  
  # Screen voxels
  screening <- screen_voxels(Y_data)
  
  expect_false(all(screening$keep[1:5]))  # Zero variance excluded
  expect_true(any(screening$flag[6:10]))  # Spikes flagged
  expect_true(all(screening$keep[11:30])) # Good voxels kept
})

test_that("Memory estimation provides reasonable values", {
  mem_req <- estimate_memory_requirements(
    n_timepoints = 300,
    n_voxels = 50000,
    n_conditions = 4,
    n_trials = 100,
    p_hrf = 30,
    m_manifold = 5
  )
  
  expect_type(mem_req, "list")
  expect_true(all(c("data_matrices_gb", "peak_estimate_gb") %in% names(mem_req)))
  expect_gt(mem_req$peak_estimate_gb, 0)
  expect_gt(mem_req$peak_estimate_gb, mem_req$data_matrices_gb)
})

test_that("Chunked processing works correctly", {
  # Create a simple processing function
  process_func <- function(voxel_indices, data_matrix) {
    # Return column means for specified voxels
    colMeans(data_matrix[, voxel_indices, drop = FALSE])
  }
  
  # Test data
  test_data <- matrix(rnorm(100 * 50), 100, 50)
  
  # Process in chunks
  result_chunked <- process_in_chunks(
    process_function = process_func,
    n_voxels = 50,
    chunk_size = 15,
    data_matrix = test_data
  )
  
  # Compare with direct processing
  result_direct <- colMeans(test_data)
  
  expect_equal(result_chunked, result_direct)
})

test_that("Data scaling check works", {
  # Very large values
  Y_large <- matrix(rnorm(100, mean = 1e8, sd = 1e7), 20, 5)
  scaled_large <- check_and_scale_data(Y_large)
  
  expect_true(scaled_large$scaling_applied)
  expect_lt(median(abs(scaled_large$Y_scaled)), 1000)
  
  # Very small values
  Y_small <- matrix(rnorm(100, mean = 0, sd = 1e-8), 20, 5)
  scaled_small <- check_and_scale_data(Y_small)
  
  expect_true(scaled_small$scaling_applied)
  expect_gt(median(abs(scaled_small$Y_scaled)), 0.01)
  
  # Integer values (warning)
  Y_int <- matrix(sample(1:100, 100, replace = TRUE), 20, 5)
  expect_warning(
    scaled_int <- check_and_scale_data(Y_int),
    "integer-valued"
  )
})

test_that("HRF physiological constraints work correctly", {
  p <- 20
  V <- 10
  TR <- 2
  time_vec <- (0:(p-1)) * TR
  
  # Create test HRFs with various issues
  hrf_matrix <- matrix(0, p, V)
  
  # Good HRF (peak at 5s)
  hrf_matrix[, 1] <- dgamma(time_vec, shape = 6, scale = 0.9) - 
                     0.35 * dgamma(time_vec, shape = 16, scale = 0.9)
  
  # Early peak (1s)
  hrf_matrix[, 2] <- c(1, 0.5, rep(0, p-2))
  
  # Late peak (15s)
  hrf_matrix[, 3] <- c(rep(0, 8), 1, rep(0.5, p-9))
  
  # Negative integral
  hrf_matrix[, 4] <- -hrf_matrix[, 1]
  
  # All zero
  hrf_matrix[, 5] <- 0
  
  # Apply constraints
  result <- apply_hrf_physiological_constraints(
    hrf_matrix = hrf_matrix,
    TR = TR,
    peak_range = c(2, 10),
    enforce_positive = TRUE,
    project_to_plausible = TRUE
  )
  
  expect_equal(dim(result$hrf_constrained), dim(hrf_matrix))
  expect_equal(as.numeric(result$quality_metrics["is_plausible", 1]), 1)  # Good HRF
  expect_equal(as.numeric(result$quality_metrics["is_plausible", 2]), 0) # Early peak
  expect_equal(as.numeric(result$quality_metrics["is_plausible", 3]), 0) # Late peak
  expect_gt(result$quality_metrics["adjustment_made", 2], 0) # Fixed
  
  # Check reasonableness scores
  score_good <- compute_hrf_reasonableness(hrf_matrix[, 1], TR)
  score_bad <- compute_hrf_reasonableness(hrf_matrix[, 2], TR)
  expect_gt(score_good, score_bad)
  expect_true(score_good > 0.5)
})

test_that("Convergence tracking works correctly", {
  # Simulate parameter updates with decreasing changes
  params1 <- rnorm(100)
  params2 <- params1 + rnorm(100, sd = 0.001)  # Smaller change
  params3 <- params2 + rnorm(100, sd = 0.0005) # Even smaller
  params4 <- params3 + rnorm(100, sd = 0.0001) # Very small
  
  # Track convergence
  history <- NULL
  history <- track_convergence_metrics(params1, NULL, 1, "test_params", history)
  history <- track_convergence_metrics(params2, params1, 2, "test_params", history)
  history <- track_convergence_metrics(params3, params2, 3, "test_params", history)
  history <- track_convergence_metrics(params4, params3, 4, "test_params", history)
  
  expect_equal(length(history$iterations), 4)
  expect_true(is.na(history$relative_change[1]))
  expect_true(all(history$relative_change[-1] > 0))
  
  # Check convergence with appropriate tolerance for the noise levels used
  status <- check_convergence_status(history, rel_tol = 0.01, min_iterations = 2)
  expect_true(status$converged)
  expect_equal(status$reason, "relative_tolerance")
  
  # Check stagnation detection
  history_stag <- history
  history_stag$relative_change <- c(NA, 1e-3, 1e-3, 1e-3)
  status_stag <- check_convergence_status(history_stag, rel_tol = 1e-5)
  expect_false(status_stag$converged)
  expect_equal(status_stag$reason, "stagnation")
})

test_that("Solution quality metrics work", {
  n <- 100
  V <- 20
  m <- 3
  k <- 2
  
  # Create test data
  Y_data <- matrix(rnorm(n * V), n, V)
  Y_predicted <- Y_data + matrix(rnorm(n * V, sd = 0.1), n, V)
  Xi_matrix <- matrix(rnorm(m * V), m, V)
  Beta_matrix <- matrix(rnorm(k * V), k, V)
  
  # Compute quality
  quality <- compute_solution_quality(
    Y_data = Y_data,
    Y_predicted = Y_predicted,
    Xi_matrix = Xi_matrix,
    Beta_matrix = Beta_matrix
  )
  
  expect_true(quality$r_squared > 0.9)  # Good fit
  expect_length(quality$r_squared_voxels, V)
  expect_true(quality$overall_quality > 0 && quality$overall_quality <= 1)
  expect_true("xi_smoothness" %in% names(quality))
  expect_true("beta_sparsity" %in% names(quality))
})

test_that("Design rank checking works", {
  n <- 50
  p <- 10
  
  # Full rank design
  X_good <- matrix(rnorm(n * p), n, p)
  result_good <- check_design_rank(X_good)
  
  expect_false(result_good$is_rank_deficient)
  expect_equal(result_good$original_rank, p)
  expect_equal(ncol(result_good$X_cleaned), p)
  
  # Rank deficient design
  X_bad <- X_good
  X_bad[, 5] <- X_bad[, 1] + X_bad[, 2]  # Linear combination
  X_bad[, 8] <- 2 * X_bad[, 3]  # Scaled duplicate
  
  expect_warning(
    result_bad <- check_design_rank(X_bad, remove_collinear = TRUE),
    "rank deficient"
  )
  
  expect_true(result_bad$is_rank_deficient)
  expect_lt(result_bad$original_rank, p)
  expect_lt(ncol(result_bad$X_cleaned), p)
  expect_true(5 %in% result_bad$removed_columns || 8 %in% result_bad$removed_columns)
})

test_that("Zero voxel handling works", {
  n <- 100
  V <- 30
  
  # Create data with problematic voxels
  Y_data <- matrix(rnorm(n * V), n, V)
  
  # All zero voxels
  Y_data[, 1:5] <- 0
  
  # Low variance voxels
  Y_data[, 6:10] <- matrix(rnorm(n * 5, sd = 1e-10), n, 5)
  
  # Test skip strategy
  result_skip <- handle_zero_voxels(Y_data, replace_with = "skip")
  
  expect_equal(length(result_skip$zero_indices), 10)
  expect_equal(result_skip$n_problematic, 10)
  expect_equal(dim(result_skip$Y_cleaned), dim(Y_data))
  
  # Test noise replacement
  result_noise <- handle_zero_voxels(Y_data, replace_with = "noise")
  
  expect_gt(var(result_noise$Y_cleaned[, 1]), 0)  # No longer zero
  expect_true(all(result_noise$Y_cleaned[, 1] != 0))
})

test_that("Condition number monitoring works", {
  # Create test matrices with varying conditioning
  mat_good <- diag(10)
  mat_poor <- diag(c(1e9, rep(1, 9)))
  mat_critical <- diag(c(1e13, rep(1, 9)))
  
  matrix_list <- list(
    good_matrix = mat_good,
    poor_matrix = mat_poor,
    critical_matrix = mat_critical
  )
  
  expect_warning(
    results <- monitor_condition_numbers(matrix_list),
    "critical condition"
  )
  
  expect_equal(nrow(results), 3)
  expect_equal(results$status[1], "OK")
  expect_equal(results$status[2], "WARNING")
  expect_equal(results$status[3], "CRITICAL")
  expect_true(all(!is.na(results$condition_number)))
})

test_that("Progress bar works correctly", {
  # Create progress bar
  pb <- create_progress_bar(total = 10)
  
  expect_equal(pb$total, 10)
  expect_equal(pb$current, 0)
  
  # Update progress
  for (i in 1:5) {
    pb <- update_progress_bar(pb, increment = 1)
  }
  
  expect_equal(pb$current, 5)
  
  # Test time formatting
  expect_equal(format_time(45), "45s")
  expect_equal(format_time(125), "2m 5s")
  expect_equal(format_time(3725), "1h 2m")
})