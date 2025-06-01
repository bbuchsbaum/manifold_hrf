library(testthat)

context("Design matrix truncation")

test_that(".create_design_matrices detects truncated HRFs", {
  events <- data.frame(
    onset = c(10, 95),
    condition = c("A", "A"),
    duration = 0
  )

  res <- manifoldhrf:::`.create_design_matrices`(
    events = events,
    n_timepoints = 100,
    TR = 1,
    hrf_length = 20
  )

  expect_equal(res$n_truncated_hrfs, 1)
  expect_equal(nrow(res$X_trial_list[[2]]), 100)
  expect_equal(ncol(res$X_trial_list[[2]]), 20)
})
