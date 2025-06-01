test_that(".determine_n_jobs handles invalid inputs", {
  expect_equal(manifoldhrf:::`.determine_n_jobs`(NULL), 1L)
  expect_equal(manifoldhrf:::`.determine_n_jobs`(0), 1L)
  expect_equal(manifoldhrf:::`.determine_n_jobs`(-5), 1L)
})
