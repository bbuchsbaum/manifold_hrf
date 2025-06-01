context("simulate_mhrf_dataset")

test_that("simulate_mhrf_dataset creates valid dataset", {
  set.seed(123)
  sim <- simulate_mhrf_dataset(
    n_voxels = 20,
    n_timepoints = 80,
    n_trials = 5,
    n_conditions = 2,
    TR = 1.5,
    hrf_variability = "moderate",
    noise_level = 5
  )

  expect_type(sim, "list")
  expect_true("dataset" %in% names(sim))
  expect_true(nrow(sim$bold_data) == 80)
  expect_true(ncol(sim$bold_data) == 20)
  expect_s3_class(sim$dataset, "matrix_dataset")
})

