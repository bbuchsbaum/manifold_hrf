context("parallel utilities")

fun <- function(x) x * x

 test_that(".parallel_lapply works sequentially and in parallel", {
  X <- 1:5
  res_seq <- manifoldhrf:::`.parallel_lapply`(X, fun, n_jobs = 1)
  res_par <- manifoldhrf:::`.parallel_lapply`(X, fun, n_jobs = 2)
  expect_equal(res_seq, res_par)
})

 test_that(".parallel_lapply handles single core gracefully", {
  X <- 1:3
  res0 <- manifoldhrf:::`.parallel_lapply`(X, fun, n_jobs = 1)
  res1 <- manifoldhrf:::`.parallel_lapply`(X, fun, n_jobs = 1)
  expect_equal(res0, res1)
})
