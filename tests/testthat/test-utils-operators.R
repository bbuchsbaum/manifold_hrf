library(testthat)

test_that("null-coalesce operator %||% works correctly", {
  expect_equal(NULL %||% "fallback", "fallback")
  expect_equal(0 %||% "fallback", 0)
  expect_equal("value" %||% "fallback", "value")
})