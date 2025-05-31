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
  L_dup[, 2] <- L_dup[, 1] * 1.001  # Near duplicate
  
  quality_dup <- check_hrf_library_quality(L_dup)
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
  expect_equal(diag(BtB), rep(1, pca_result$m_final_dim), tolerance = 1e-10)
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
  
  # Test aggressive preset
  params_agg <- get_preset_params("aggressive")
  expect_gt(params_agg$m_manifold_dim_target, 5)
  expect_lt(params_agg$lambda_gamma, 0.01)
  
  # Test with data scaling
  params_scaled <- get_preset_params("balanced", data_scale = 10)
  expect_gt(params_scaled$lambda_gamma, get_preset_params("balanced")$lambda_gamma)
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
  result <- extract_xi_beta_raw_svd_robust(
    Gamma_coeffs_matrix = Gamma,
    m_manifold_dim = m,
    k_conditions = k
  )
  
  expect_equal(dim(result$Xi_raw_matrix), c(m, V))
  expect_equal(dim(result$Beta_raw_matrix), c(k, V))
  expect_false(any(is.na(result$Xi_raw_matrix)))
  expect_false(any(is.infinite(result$Xi_raw_matrix)))
  
  # Check quality metrics
  expect_equal(result$quality_metrics$svd_method[1], "zero")
  expect_true(any(result$quality_metrics$regularization_applied))
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