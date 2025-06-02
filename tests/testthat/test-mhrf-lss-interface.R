# Tests for mhrf_lss user interface

test_that("mhrf_lss interface validates inputs correctly", {
  # Create minimal test data
  set.seed(123)
  n <- 100
  V <- 20
  
  # Mock dataset structure
  mock_dataset <- list(
    TR = 2,
    run_length = c(50, 50),
    event_table = data.frame(
      onset = c(10, 30, 60, 80),
      duration = 1,
      condition = c("A", "B", "A", "B"),
      run = c(1, 1, 2, 2),
      block = c(1, 1, 2, 2)
    )
  )
  
  # Mock get_data_matrix function
  get_data_matrix <- function(dataset) {
    matrix(rnorm(100 * 20), 100, 20)
  }
  
  # Test that it requires fmrireg
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    expect_error(
      manifoldhrf:::mhrf_lss(~ hrf(condition), dataset = mock_dataset),
      "fmrireg"
    )
  }
})

test_that("create_hrf_manifold handles different input types", {
  skip_if_not_installed("fmrireg")
  # Test with preset parameters
  custom_params <- get_preset_params("balanced")
  custom_params$k_local_nn_for_sigma <- 2  # Small value for canonical (only 3 HRFs)
  expect_silent(
    manifold1 <- suppressMessages(suppressWarnings(create_hrf_manifold(
      hrf_library = "canonical",
      params = custom_params,
      TR = 2,
      verbose = FALSE
    )))
  )
  
  # Test with custom parameters
  custom_params2 <- list(
    m_manifold_dim_target = 5,
    k_local_nn_for_sigma = 2,  # Small value for canonical (only 3 HRFs)
    TR_precision = 0.1
  )
  
  expect_silent(
    manifold2 <- suppressMessages(suppressWarnings(create_hrf_manifold(
      hrf_library = "canonical",
      params = custom_params2,
      TR = 2,
      verbose = FALSE
    )))
  )
  
  # Test with matrix input
  hrf_matrix <- matrix(rnorm(20 * 10), 20, 10)
  expect_silent(
    manifold3 <- suppressMessages(suppressWarnings(create_hrf_manifold(
      hrf_library = hrf_matrix,
      params = "balanced",
      TR = 2,
      verbose = FALSE
    )))
  )
})

test_that("preset parameter retrieval works", {
  # Test all presets
  presets <- c("conservative", "balanced", "aggressive")
  
  for (preset in presets) {
    params <- get_preset_params(preset)
    expect_type(params, "list")
    expect_true("m_manifold_dim_target" %in% names(params))
    expect_true("lambda_gamma" %in% names(params))
  }
  
  # Test with data scaling
  params_scaled <- get_preset_params("balanced", data_scale = 10)
  params_base <- get_preset_params("balanced", data_scale = 1)
  
  expect_gt(params_scaled$lambda_gamma, params_base$lambda_gamma)
})

test_that("create_trial_matrices generates correct dimensions", {
  # Skip this test if fmrireg is not available
  skip_if_not_installed("fmrireg")
  
  # This test would require actual fmrireg objects
  # For now, we skip it as it requires deep integration with fmrireg
  skip("Requires fmrireg integration")
})

test_that("run_mhrf_lss_chunked handles chunking correctly", {
  set.seed(456)
  n <- 100
  V <- 50
  k <- 3
  m <- 4
  p <- 20
  
  # Create test data
  Y_data <- matrix(rnorm(n * V), n, V)
  
  # Mock design info
  design_info <- list(
    X_condition_list = lapply(1:k, function(i) matrix(rnorm(n * p), n, p)),
    X_trial_list = lapply(1:10, function(i) matrix(rnorm(n * p), n, p)),
    n_conditions = k,
    n_trials = 10
  )
  
  # Mock manifold
  manifold <- list(
    B_reconstructor_matrix = matrix(rnorm(p * m), p, m),
    m_manifold_dim = m,
    parameters = get_preset_params("balanced")
  )
  
  # Mock confounds
  Z_confounds <- matrix(rnorm(n * 2), n, 2)
  
  # Test with different chunk sizes
  for (nchunks in c(1, 5, 10)) {
    result <- run_mhrf_lss_chunked(
      Y_data = Y_data,
      design_info = design_info,
      manifold = manifold,
      Z_confounds = Z_confounds,
      voxel_coords = NULL,
      params = manifold$parameters,
      outlier_weights = NULL,
      estimation = "condition",
      nchunks = nchunks,
      progress = FALSE
    )
    
    expect_equal(ncol(result$H_shapes), V)
    expect_equal(ncol(result$Xi_smoothed), V)
    expect_equal(ncol(result$Beta_condition), V)
    expect_equal(result$diagnostics$n_chunks, nchunks)
  }
})

test_that("S3 methods work for mhrf_lss_result", {
  # Create mock fit object
  mock_fit <- structure(
    list(
      parameters = list(manifold = list(), robust = FALSE),
      manifold = list(
        method_used = "diffusion_maps",
        m_manifold_dim = 5,
        parameters = list(n_hrfs_library = 30)
      ),
      hrf = list(
        raw = matrix(rnorm(20 * 100), 20, 100),
        smoothed = matrix(rnorm(20 * 100), 20, 100)
      ),
      beta = list(
        condition_initial = matrix(rnorm(3 * 100), 3, 100),
        condition_final = matrix(rnorm(3 * 100), 3, 100),
        trial = matrix(rnorm(50 * 100), 50, 100)
      ),
      qc = list(),
      diagnostics = list(),
      call = quote(mhrf_lss())
    ),
    class = "mhrf_lss_result"
  )
  
  # Test print method
  expect_output(print(mock_fit), "M-HRF-LSS Result")
  
  # Test coef method
  betas_cond <- coef(mock_fit, type = "condition")
  expect_equal(dim(betas_cond), c(3, 100))
  
  betas_trial <- coef(mock_fit, type = "trial")
  expect_equal(dim(betas_trial), c(50, 100))
  
  hrfs <- coef(mock_fit, type = "hrf")
  expect_equal(dim(hrfs), c(20, 100))
  
  # Test plot method (just check it doesn't error)
  expect_error(
    plot(mock_fit, type = "hrfs", voxels = 1:3),
    NA
  )
})

test_that("extract_voxel_coordinates handles different input types", {
  # Test with matrix coordinates
  coords_matrix <- matrix(1:30, 10, 3)
  result1 <- manifoldhrf:::extract_voxel_coordinates(coords_matrix, mask = NULL)
  expect_equal(result1, coords_matrix)
  
  # Test with NULL
  result2 <- manifoldhrf:::extract_voxel_coordinates(NULL, mask = NULL)
  expect_null(result2)
  
  # Test with mask array
  mask <- array(0, dim = c(5, 5, 2))
  mask[1:3, 1:3, 1] <- 1
  result3 <- manifoldhrf:::extract_voxel_coordinates(mask, mask = mask)
  expect_equal(ncol(result3), 3)
  expect_equal(nrow(result3), sum(mask != 0))
})

test_that("compute_pipeline_diagnostics returns expected structure", {
  Y_data <- matrix(rnorm(100 * 50), 100, 50)
  Y_proj <- Y_data
  H_shapes <- matrix(rnorm(20 * 50), 20, 50)
  Beta_condition <- matrix(rnorm(3 * 50), 3, 50)
  Xi_smoothed <- matrix(rnorm(5 * 50), 5, 50)
  
  diag <- compute_pipeline_diagnostics(
    Y_data, Y_proj, H_shapes, Beta_condition, Xi_smoothed
  )
  
  expect_type(diag, "list")
  expect_true("n_voxels_processed" %in% names(diag))
  expect_true("manifold_variance" %in% names(diag))
  expect_true("hrf_stats" %in% names(diag))
  expect_true("timestamp" %in% names(diag))
  
  expect_equal(diag$n_voxels_processed, 50)
  expect_length(diag$manifold_variance, 5)
})