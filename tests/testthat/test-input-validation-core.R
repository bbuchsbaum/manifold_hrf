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
# Priority 1 tests from TESTING_PRIORITIES.md

 test_that(".validate_events errors when events exceed data length", {
  events <- data.frame(
    onset = c(0, 20),
    duration = 1,
    condition = c("A", "B")
  )
  expect_error(
    manifoldhrf:::`.validate_events`(events, n_timepoints = 5, TR = 2),
    "after data ends"
  )
 })

 test_that(".validate_Y_data errors on non-finite values", {
  Y <- matrix(c(rep(NA, 18), rep(1, 12)), nrow = 15, ncol = 2)
  expect_error(
    manifoldhrf:::`.validate_Y_data`(Y),
    "non-finite values"
  )
 })

 test_that(".validate_events preserves inconsistent factor levels", {
  events <- data.frame(
    onset = c(0, 5, 10),
    duration = 1,
    condition = c("face", "Face", "face ")
  )
  info <- manifoldhrf:::`.validate_events`(events, n_timepoints = 20, TR = 1)
  expect_equal(info$n_conditions, 3)
  expect_setequal(info$conditions, c("face", "Face", "face "))
 })
