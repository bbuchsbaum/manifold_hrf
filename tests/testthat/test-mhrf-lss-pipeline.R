context("mhrf_lss pipeline integration")

library(testthat)


test_that("LSS failure propagation does not crash pipeline", {
  skip_if_not_installed("fmrilss")

  set.seed(1)
  n <- 20
  p <- 6
  V <- 3

  Y <- matrix(rnorm(n * V), n, V)

  # Create two identical trial regressors to induce rank deficiency
  X_base <- matrix(0, n, p)
  if (n >= p + 4) {
    X_base[5:(5 + p - 1), ] <- diag(p)
  }
  X_trials <- list(X_base, X_base)

  design_info <- list(
    X_condition_list = list(matrix(rnorm(n * p), n, p)),
    X_trial_list = X_trials,
    n_conditions = 1,
    n_trials = length(X_trials)
  )

  params <- get_preset_params("balanced")
  params$lambda_spatial_smooth <- 0
  manifold <- list(
    B_reconstructor_matrix = matrix(rnorm(p * 2), p, 2),
    m_manifold_dim = 2,
    library_hrfs = matrix(1, p, 1),
    parameters = params
  )

  Z_confounds <- matrix(0, n, 1)

  expect_error(
    fit <- run_mhrf_lss_standard(
      Y_data = Y,
      design_info = design_info,
      manifold = manifold,
      Z_confounds = Z_confounds,
      voxel_coords = NULL,
      params = manifold$parameters,
      outlier_weights = NULL,
      estimation = "trial",
      progress = FALSE
    ),
    NA
  )

  expect_equal(dim(fit$Beta_trial), c(length(X_trials), V))
})


test_that("inconsistent masking between stages triggers error", {
  skip_if_not_installed("fmrilss")

  set.seed(2)
  n <- 20
  p <- 6
  V <- 5

  Y <- matrix(rnorm(n * V), n, V)

  design_info <- list(
    X_condition_list = list(matrix(rnorm(n * p), n, p)),
    X_trial_list = list(matrix(rnorm(n * p), n, p)),
    n_conditions = 1,
    n_trials = 1
  )

  params <- get_preset_params("balanced")
  params$num_neighbors_Lsp <- 2
  manifold <- list(
    B_reconstructor_matrix = matrix(rnorm(p * 2), p, 2),
    m_manifold_dim = 2,
    library_hrfs = matrix(1, p, 1),
    parameters = params
  )

  Z_confounds <- matrix(0, n, 1)

  # Provide coordinates for fewer voxels than present in the data
  voxel_coords <- matrix(runif((V - 1) * 3), ncol = 3)

  expect_error(
    run_mhrf_lss_standard(
      Y_data = Y,
      design_info = design_info,
      manifold = manifold,
      Z_confounds = Z_confounds,
      voxel_coords = voxel_coords,
      params = manifold$parameters,
      outlier_weights = NULL,
      estimation = "condition",
      progress = FALSE
    )
  )
})


test_that("metadata alignment handles non sequential trials", {
  skip_if_not_installed("fmrilss")

  set.seed(3)
  n <- 20
  p <- 5
  V <- 4

  Y <- matrix(rnorm(n * V), n, V)

  X_trials <- replicate(3, matrix(rnorm(n * p), n, p), simplify = FALSE)
  event_table <- data.frame(
    onset = seq(0, by = 2, length.out = 4)[-3],
    duration = 0,
    condition = "A",
    trial_id = c(1, 2, 4)
  )

  design_info <- list(
    X_condition_list = list(matrix(rnorm(n * p), n, p)),
    X_trial_list = X_trials,
    event_table = event_table,
    n_conditions = 1,
    n_trials = length(X_trials)
  )

  params <- get_preset_params("balanced")
  params$lambda_spatial_smooth <- 0
  manifold <- list(
    B_reconstructor_matrix = matrix(rnorm(p * 2), p, 2),
    m_manifold_dim = 2,
    library_hrfs = matrix(1, p, 1),
    parameters = params
  )

  Z_confounds <- matrix(0, n, 1)

  fit <- run_mhrf_lss_standard(
    Y_data = Y,
    design_info = design_info,
    manifold = manifold,
    Z_confounds = Z_confounds,
    voxel_coords = NULL,
    params = manifold$parameters,
    outlier_weights = NULL,
    estimation = "trial",
    progress = FALSE
  )

  expect_equal(nrow(fit$Beta_trial), nrow(event_table))
})

