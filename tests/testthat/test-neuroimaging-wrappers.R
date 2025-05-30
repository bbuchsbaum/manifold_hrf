# Tests for Neuroimaging Wrapper Functions (EPIC 5)
# Tests for MHRF-NIM-* interface expectations

library(testthat)

# The wrapper functions are not yet implemented, so they should error
# with the message "Function not yet implemented" when called.

test_that("construct_hrf_manifold_nim has correct interface", {
  expect_error(
    construct_hrf_manifold_nim(NULL, 0.5),
    "Function not yet implemented"
  )
})

test_that("process_subject_mhrf_lss_nim and package_mhrf_results_nim have correct interface", {
  expect_error(
    process_subject_mhrf_lss_nim(NULL, NULL, NULL, NULL, NULL, list()),
    "Function not yet implemented"
  )
  expect_error(
    package_mhrf_results_nim(list(), NULL, NULL, list(), list()),
    "Function not yet implemented"
  )
})
