# Tests for Neuroimaging Wrapper Functions (EPIC 5)
# Tests for MHRF-NIM-* interface expectations

library(testthat)

test_that("construct_hrf_manifold_nim returns basis list", {
  L <- matrix(seq(0, 4, length.out = 6), nrow = 2)
  res <- construct_hrf_manifold_nim(L, 0.5)
  expect_type(res, "list")
  expect_true("B_reconstructor_matrix" %in% names(res))
})

test_that("process_subject_mhrf_lss_nim and package_mhrf_results_nim return lists", {
  dummy_bold <- matrix(rnorm(20), nrow = 5)
  manifold <- construct_hrf_manifold_nim(matrix(seq(0, 4, length.out = 6), nrow = 2), 0.5)
  res_subj <- process_subject_mhrf_lss_nim(dummy_bold, NULL, NULL, NULL, manifold, list())
  expect_type(res_subj, "list")
  res_out <- package_mhrf_results_nim(list(), NULL, NULL, list(), list())
  expect_s3_class(res_out, "mhrf_results")
})
