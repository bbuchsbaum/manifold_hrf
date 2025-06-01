context("Output directory handling")

# Test that output directory is created when missing

test_that("mhrf_analyze creates missing output directory", {
  set.seed(1)
  n <- 20
  V <- 5
  Y_data <- matrix(rnorm(n * V), n, V)

  events <- data.frame(
    condition = "A",
    onset = c(5, 10),
    duration = 1
  )

  out_dir <- file.path(tempdir(), "mhrf_output_test")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  expect_false(dir.exists(out_dir))

  result <- mhrf_analyze(
    Y_data = Y_data,
    events = events,
    TR = 2,
    save_intermediate = TRUE,
    output_dir = out_dir,
    verbose = 0
  )

  expect_s3_class(result, "mhrf_result")
  expect_true(dir.exists(out_dir))
  expect_true(file.exists(file.path(out_dir, "mhrf_result.rds")))
})
