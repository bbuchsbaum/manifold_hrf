context("input validation core")

 test_that(".validate_Y_data rejects invalid input", {
  expect_error(manifoldhrf:::`.validate_Y_data`("bad"), "matrix or NeuroVec")
 })

 test_that(".validate_events detects missing columns", {
  events <- data.frame(onset = 1:5)
  expect_error(manifoldhrf:::`.validate_events`(events, n_timepoints = 10, TR = 1),
               "missing required columns")
 })

 test_that(".validate_voxel_mask checks length", {
  mask <- c(TRUE, FALSE)
  expect_error(manifoldhrf:::`.validate_voxel_mask`(mask, n_voxels = 5),
               "length")
 })

 test_that(".validate_parameters validates TR and preset", {
  expect_error(manifoldhrf:::`.validate_parameters`(TR = -1, preset = "balanced", n_voxels = 10),
               "TR value seems unrealistic")
  expect_error(manifoldhrf:::`.validate_parameters`(TR = 2, preset = "unknown", n_voxels = 10),
               "Invalid preset")
 })
