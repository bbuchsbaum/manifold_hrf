library(testthat)

test_that("select_manifold_dim chooses correct dimension", {
  eig <- c(1, 0.5, 0.3, 0.1, 0.05)
  res <- select_manifold_dim(eig, min_var = 0.80)
  expect_equal(res$m_auto, 2)
  expect_true(res$cum_var[res$m_auto] >= 0.80)
})

test_that("select_manifold_dim handles edge cases", {
  expect_warning(res <- select_manifold_dim(c(1, 0, 0), 0.9))
  expect_equal(res$m_auto, 1)
  expect_error(select_manifold_dim(c(1, 0.1), 1.5), "between 0 and 1")
})
