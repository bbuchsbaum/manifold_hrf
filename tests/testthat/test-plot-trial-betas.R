library(testthat)

# Simple object with trial amplitudes
mock_result <- structure(
  list(trial_amplitudes = matrix(rnorm(10), 5, 2)),
  class = "mhrf_result"
)


test_that("plot_trial_betas runs without error", {
  expect_silent(plot_trial_betas(mock_result, voxel = 1))
})
