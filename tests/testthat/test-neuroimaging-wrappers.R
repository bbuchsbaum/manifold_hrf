# Tests for Neuroimaging Wrapper Functions
# Tests for MHRF-NIM-IO-MANIFOLD-01

test_that("construct_hrf_manifold_nim works with predefined libraries", {
  # Test FLOBS library - expect warning due to small library size
  expect_warning(
    manifold_flobs <- construct_hrf_manifold_nim(
      hrf_library_source = "FLOBS",
      TR_precision = 0.5,
      hrf_duration = 20,
      m_manifold_dim_target = 3
    ),
    "Adjusting k_local_nn_for_sigma"
  )
  
  # Check output structure
  expect_s3_class(manifold_flobs, "mhrf_manifold")
  expect_type(manifold_flobs, "list")
  
  # Check required components
  expect_true("B_reconstructor_matrix" %in% names(manifold_flobs))
  expect_true("manifold_hrf_basis" %in% names(manifold_flobs))
  expect_true("Phi_coords_matrix" %in% names(manifold_flobs))
  expect_true("eigenvalues_S" %in% names(manifold_flobs))
  
  # Check dimensions
  expect_equal(ncol(manifold_flobs$B_reconstructor_matrix), 
               manifold_flobs$m_manifold_dim)
  expect_equal(nrow(manifold_flobs$B_reconstructor_matrix),
               length(seq(0, 20, by = 0.5)))
  
  # Check library info
  expect_equal(manifold_flobs$library_info$name, "FLOBS")
})

test_that("construct_hrf_manifold_nim works with half-cosine library", {
  # Test half-cosine library
  manifold_hc <- construct_hrf_manifold_nim(
    hrf_library_source = "half_cosine",
    TR_precision = 0.2,
    m_manifold_dim_target = 4
  )
  
  expect_s3_class(manifold_hc, "mhrf_manifold")
  expect_equal(manifold_hc$library_info$name, "half_cosine")
  expect_lte(manifold_hc$m_manifold_dim, 4)
})

test_that("construct_hrf_manifold_nim works with gamma grid library", {
  # Test gamma grid library
  manifold_gg <- construct_hrf_manifold_nim(
    hrf_library_source = "gamma_grid",
    TR_precision = 0.1,
    m_manifold_dim_target = 5
  )
  
  expect_s3_class(manifold_gg, "mhrf_manifold")
  expect_equal(manifold_gg$library_info$name, "gamma_grid")
  expect_equal(manifold_gg$parameters$n_hrfs_library, 100)  # 10x10 grid
})

test_that("construct_hrf_manifold_nim handles custom matrix input", {
  # Create custom HRF matrix
  time_points <- seq(0, 24, by = 0.5)
  p <- length(time_points)
  N <- 20
  
  # Create some HRFs
  L_custom <- matrix(0, p, N)
  for (i in 1:N) {
    L_custom[, i] <- dgamma(time_points, shape = 4 + i/5, rate = 1)
  }
  L_custom <- apply(L_custom, 2, function(x) x / sum(x))
  
  manifold_custom <- construct_hrf_manifold_nim(
    hrf_library_source = L_custom,
    TR_precision = 0.5,
    m_manifold_dim_target = 4
  )
  
  expect_s3_class(manifold_custom, "mhrf_manifold")
  expect_equal(manifold_custom$library_info$name, "custom_matrix")
  expect_equal(manifold_custom$parameters$n_hrfs_library, N)
})

test_that("construct_hrf_manifold_nim validates inputs", {
  # Invalid library source
  expect_error(
    construct_hrf_manifold_nim("invalid_library"),
    "Unknown HRF library source"
  )
  
  # File that doesn't exist
  expect_error(
    construct_hrf_manifold_nim("/path/to/nonexistent/file.rds"),
    "file not found"
  )
  
  # Invalid matrix dimensions
  bad_matrix <- matrix(1:5, 1, 5)  # Only 1 time point
  expect_error(
    construct_hrf_manifold_nim(bad_matrix),
    "must have at least 2 time points"
  )
})

test_that("construct_hrf_manifold_nim handles sparse matrices for large libraries", {
  # Create large library
  time_points <- seq(0, 20, by = 0.5)
  p <- length(time_points)
  N <- 6000  # Above sparse threshold
  
  # Create HRFs (simplified for speed)
  L_large <- matrix(rnorm(p * N), p, N)
  L_large <- apply(L_large, 2, function(x) x / sum(abs(x)))
  
  # Should use sparse matrices
  expect_message(
    manifold_large <- construct_hrf_manifold_nim(
      hrf_library_source = L_large,
      TR_precision = 0.5,
      m_manifold_dim_target = 5,
      sparse_threshold = 5000
    ),
    "Using sparse matrices"
  )
  
  expect_true(manifold_large$parameters$use_sparse)
})

test_that("construct_hrf_manifold_nim adjusts target dimension when needed", {
  # Create small library
  time_points <- seq(0, 20, by = 1)
  p <- length(time_points)
  N <- 4  # Smaller than target dimension
  
  L_small <- matrix(rnorm(p * N), p, N)
  
  # Should warn and adjust
  expect_warning(
    manifold_small <- construct_hrf_manifold_nim(
      hrf_library_source = L_small,
      m_manifold_dim_target = 5  # Larger than N-1
    ),
    "Adjusting target"
  )
  
  expect_lte(manifold_small$m_manifold_dim, N - 1)
})

test_that("print and summary methods work for mhrf_manifold", {
  manifold <- construct_hrf_manifold_nim(
    hrf_library_source = "gamma_grid",
    TR_precision = 0.5,
    m_manifold_dim_target = 4
  )
  
  # Test print method
  expect_output(print(manifold), "M-HRF Manifold Object")
  expect_output(print(manifold), "Library: gamma_grid")
  expect_output(print(manifold), "Variance explained:")
  
  # Test summary method
  expect_output(summary(manifold), "Eigenvalue spectrum")
  expect_output(summary(manifold), "Reconstructor matrix dimensions:")
})

test_that("manifold HRF basis object is created correctly", {
  expect_warning(
    manifold <- construct_hrf_manifold_nim(
      hrf_library_source = "FLOBS",
      TR_precision = 0.5,
      m_manifold_dim_target = 3
    ),
    "Adjusting k_local_nn_for_sigma"
  )
  
  hrf_basis <- manifold$manifold_hrf_basis
  
  # Check structure
  expect_type(hrf_basis, "list")
  expect_s3_class(hrf_basis, "mhrf_basis")
  expect_equal(hrf_basis$type, "manifold")
  expect_equal(hrf_basis$nbasis, manifold$m_manifold_dim)
  expect_true("evaluate" %in% names(hrf_basis))
  expect_true(is.function(hrf_basis$evaluate))
})

test_that("construct_hrf_manifold_nim handles list input", {
  # Test with list of HRF objects
  # For now, this creates dummy HRFs since we don't have fmrireg loaded
  hrf_list <- list(
    hrf1 = list(name = "hrf1"),
    hrf2 = list(name = "hrf2"),
    hrf3 = list(name = "hrf3")
  )
  
  expect_warning(
    manifold_list <- construct_hrf_manifold_nim(
      hrf_library_source = hrf_list,
      TR_precision = 0.5,
      m_manifold_dim_target = 2
    ),
    "Adjusting k_local_nn_for_sigma"
  )
  
  expect_s3_class(manifold_list, "mhrf_manifold")
  expect_equal(manifold_list$library_info$name, "custom_list")
  expect_equal(manifold_list$library_info$n_hrfs, 3)
})

test_that("placeholder functions raise appropriate errors", {
  # Test that unimplemented functions give clear errors
  expect_error(
    process_subject_mhrf_lss_nim(NULL, NULL, NULL),
    "not yet implemented"
  )
  
  expect_error(
    package_mhrf_results_nim(NULL, NULL, NULL, NULL, NULL),
    "not yet implemented"
  )
})
