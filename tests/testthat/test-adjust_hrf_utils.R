library(testthat)

# Tests for adjust_hrf_for_bounds and .validate_and_standardize_lambda


test_that("adjust_hrf_for_bounds truncates and pads correctly", {
  hrf <- 1:5
  expect_warning(trunc <- adjust_hrf_for_bounds(hrf, 3), "HRF truncated")
  expect_equal(trunc, 1:3)

  padded <- adjust_hrf_for_bounds(hrf, 7)
  expect_equal(length(padded), 7)
  expect_equal(padded[1:5], hrf)
  expect_true(all(padded[6:7] == 0))

  expect_error(adjust_hrf_for_bounds("bad", 3), "numeric")
})

test_that(".validate_and_standardize_lambda validates input", {
  expect_error(manifoldhrf:::`.validate_and_standardize_lambda`(-1, "lambda"),
               "non-negative")
  expect_error(manifoldhrf:::`.validate_and_standardize_lambda`(c(1,2), "lambda"),
               "non-negative")
  expect_equal(manifoldhrf:::`.validate_and_standardize_lambda`(0.5, "lambda"),
               0.5)
  expect_equal(manifoldhrf:::`.validate_and_standardize_lambda`(2L, "lambda"), 2)
})
