context("Taylor series utilities")

# Test numeric_taylor_coefficients and evaluate_taylor

test_that("Taylor series approximates sin near zero", {
  f <- sin
  coeff <- numeric_taylor_coefficients(f, 0, order = 5)
  x <- seq(-0.1, 0.1, length.out = 5)
  approx <- evaluate_taylor(coeff, 0, x)
  expect_equal(approx, sin(x), tolerance = 1e-6)
})
