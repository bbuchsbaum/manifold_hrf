context("fmrireg helper utilities")

# Test HRF_RAW_EVENT_BASIS -----------------------------------------------------

test_that("HRF_RAW_EVENT_BASIS creates valid HRF basis", {
  skip_if_not_installed("fmrireg")

  hrf <- manifoldhrf:::HRF_RAW_EVENT_BASIS(p_length = 4, TR_sample = 1)
  expect_s3_class(hrf, "HRF")

  # Evaluate basis at sample times
  mat <- hrf(seq(0, 3, by = 1))
  expect_equal(mat, diag(4))
  expect_equal(attr(hrf, "nbasis"), 4L)
  expect_equal(attr(hrf, "span"), 4)
})

# Test as_fmrireg_hrfs ---------------------------------------------------------

test_that("as_fmrireg_hrfs converts matrices and mhrf_result objects", {
  skip_if_not_installed("fmrihrf")
  
  mat <- matrix(rnorm(20), 5, 4)
  hrfs <- manifoldhrf::as_fmrireg_hrfs(mat, TR = 2, prefix = "vox")
  expect_length(hrfs, 4)
  expect_s3_class(hrfs[[1]], "HRF")
  expect_equal(attr(hrfs[[1]], "span"), 8)

  result_obj <- structure(
    list(
      hrf_shapes = mat,
      metadata = list(parameters = list(TR = 2))
    ),
    class = "mhrf_result"
  )
  hrfs2 <- manifoldhrf::as_fmrireg_hrfs(result_obj, prefix = "vox")
  expect_equal(length(hrfs2), 4)
})

# Test create_logger -----------------------------------------------------------

test_that("create_logger records and prints messages", {
  logger <- manifoldhrf::create_logger()
  expect_s3_class(logger, "mhrf_logger")

  logger$add("first message")
  logger$add("second message")
  msgs <- logger$get()
  expect_equal(length(msgs), 2)
  expect_true(all(grepl("message", msgs)))
  expect_output(print(logger), "M-HRF-LSS Log")
})

