library(testthat)

skip_if_not_installed("hrfals")

# simple matrices for testing
X <- matrix(1:9, 3, 3)
Y <- matrix(1:12, 3, 4)

# expect that hrfals functions return matrices of dim ncol(X) x ncol(Y)

# run_fastlss
res_a <- run_fastlss(X, Y)

test_that("run_fastlss returns matrix with correct dims", {
  expect_true(is.matrix(res_a))
  expect_equal(dim(res_a), c(ncol(X), ncol(Y)))
})

# run_stablss
res_b <- run_stablss(X, Y)

test_that("run_stablss returns matrix with correct dims", {
  expect_true(is.matrix(res_b))
  expect_equal(dim(res_b), c(ncol(X), ncol(Y)))
})
