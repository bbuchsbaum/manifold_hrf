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


test_that("simulate_mhrf_dataset respects variability and noise parameters", {
  set.seed(42)
  sim_none <- simulate_mhrf_dataset(
    n_voxels = 5,
    n_timepoints = 40,
    n_trials = 2,
    n_conditions = 2,
    hrf_variability = "none",
    noise_level = 0
  )

  hrfs <- sim_none$ground_truth$hrfs$matrix
  diffs <- apply(hrfs, 2, function(col) max(abs(col - hrfs[,1])))
  expect_true(all(diffs < 1e-8))

  dvars_zero <- manifoldhrf:::`compute_dvars`(sim_none$bold_data)

  set.seed(42)
  sim_noise <- simulate_mhrf_dataset(
    n_voxels = 5,
    n_timepoints = 40,
    n_trials = 2,
    n_conditions = 2,
    hrf_variability = "none",
    noise_level = 5
  )
  dvars_noise <- manifoldhrf:::`compute_dvars`(sim_noise$bold_data)

  expect_gt(dvars_noise, dvars_zero)
})
