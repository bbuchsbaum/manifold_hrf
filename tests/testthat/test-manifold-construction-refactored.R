# Tests for refactored manifold construction
library(testthat)
library(manifoldhrf)

test_that("build_manifold_affinity creates proper sparse affinity matrix", {
  set.seed(123)
  # Small test case
  p <- 20
  N <- 50
  L <- matrix(rnorm(p * N), p, N)
  
  # Test symmetric normalization
  aff_sym <- build_manifold_affinity(
    L, 
    k = 10,
    normalization = "symmetric",
    return_diagnostics = TRUE
  )
  
  expect_s3_class(aff_sym, "manifold_affinity")
  expect_s4_class(aff_sym$W, "sparseMatrix")
  expect_s4_class(aff_sym$T_norm, "sparseMatrix")
  expect_equal(aff_sym$normalization, "symmetric")
  
  # Check symmetry
  expect_true(Matrix::isSymmetric(aff_sym$T_norm))
  
  # Check diagnostics
  expect_true(!is.null(aff_sym$diagnostics))
  expect_equal(length(aff_sym$diagnostics$sigma_values), N)
})

test_that("build_manifold_affinity handles different normalizations", {
  set.seed(456)
  p <- 15
  N <- 30
  L <- matrix(rnorm(p * N), p, N)
  
  # Row-stochastic normalization
  aff_row <- build_manifold_affinity(
    L, 
    k = 5,
    normalization = "row_stochastic"
  )
  
  # Check row sums equal 1
  row_sums <- Matrix::rowSums(aff_row$T_norm)
  expect_true(all(abs(row_sums - 1) < 1e-10))
  
  # Symmetric normalization
  aff_sym <- build_manifold_affinity(
    L, 
    k = 5,
    normalization = "symmetric"
  )
  
  # Should be symmetric
  expect_true(Matrix::isSymmetric(aff_sym$T_norm))
  
  # Row sums should NOT equal 1 for symmetric
  row_sums_sym <- Matrix::rowSums(aff_sym$T_norm)
  expect_false(all(abs(row_sums_sym - 1) < 1e-10))
})

test_that("build_manifold_affinity handles isolated nodes correctly", {
  set.seed(789)
  # Create data with guaranteed isolated nodes
  p <- 10
  N <- 10
  
  # Main cluster of points
  L_main <- matrix(rnorm(p * (N-1)), p, N-1)
  # One point very far away
  L_isolated <- matrix(rnorm(p * 1) + 10000, p, 1)  # Extremely far
  L <- cbind(L_main, L_isolated)
  
  # With k=2, the isolated point should have no neighbors
  expect_warning(
    aff_isolate <- build_manifold_affinity(
      L, 
      k = 2,
      handle_isolated = "connect_to_self"
    ),
    "isolated nodes"
  )
  
  # Test error on isolated
  expect_error(
    build_manifold_affinity(
      L, 
      k = 2,
      handle_isolated = "error"
    ),
    "isolated nodes"
  )
  
  # Test removal
  expect_warning(
    aff_remove <- build_manifold_affinity(
      L, 
      k = 2,
      handle_isolated = "warn_remove"
    ),
    "isolated nodes"
  )
  
  # Removed matrix should be smaller
  expect_lt(nrow(aff_remove$T_norm), N)
})

test_that("build_manifold_affinity uses correct k-NN backend", {
  skip_if_not(exists("knn_search_cpp", mode = "function"))
  
  set.seed(111)
  p <- 10
  N <- 25
  L <- matrix(rnorm(p * N), p, N)
  
  # Test exact method
  aff_exact <- build_manifold_affinity(
    L, 
    k = 5,
    ann_method = "exact",
    return_diagnostics = TRUE
  )
  
  expect_equal(aff_exact$diagnostics$knn_backend, "knn_search_cpp")
  
  # Test RANN if available
  if (requireNamespace("RANN", quietly = TRUE)) {
    aff_rann <- build_manifold_affinity(
      L, 
      k = 5,
      ann_method = "RANN",
      return_diagnostics = TRUE
    )
    expect_equal(aff_rann$diagnostics$knn_backend, "RANN")
  }
})

test_that("compute_diffusion_basis works with symmetric normalization", {
  set.seed(222)
  p <- 15
  N <- 40
  L <- matrix(rnorm(p * N), p, N)
  
  # Build affinity
  aff <- build_manifold_affinity(
    L, 
    k = 8,
    normalization = "symmetric"
  )
  
  # Compute basis
  basis <- compute_diffusion_basis(
    aff,
    n_dims = 5,
    return_basis = TRUE,
    L_library_matrix = L
  )
  
  expect_s3_class(basis, "diffusion_basis")
  expect_equal(basis$n_dims, 5)
  expect_equal(dim(basis$eigenvectors), c(N, 5))
  expect_true(all(basis$eigenvalues >= 0))  # Should be non-negative for symmetric
  
  # Check reconstructor
  expect_equal(dim(basis$B_reconstructor), c(p, 5))
})

test_that("compute_diffusion_basis handles auto dimension selection", {
  set.seed(333)
  p <- 20
  N <- 60
  
  # Create data with clear low-dimensional structure
  # 3D manifold embedded in high dimensions
  t <- seq(0, 2*pi, length.out = N)
  coords_3d <- cbind(cos(t), sin(t), t/2)
  
  # Embed in high dimensions with noise
  embedding <- matrix(rnorm(p * 3), p, 3)
  L <- embedding %*% t(coords_3d) + 0.1 * matrix(rnorm(p * N), p, N)
  
  # Build affinity
  aff <- build_manifold_affinity(L, k = 10)
  
  # Auto dimension selection
  basis_auto <- compute_diffusion_basis(
    aff,
    n_dims = "auto",
    min_variance = 0.9,
    return_basis = FALSE
  )
  
  # Should select a reasonable number of dimensions
  expect_true(basis_auto$n_dims >= 2)
  expect_true(basis_auto$n_dims <= 20)
})

test_that("Integration: full pipeline with both normalizations", {
  set.seed(444)
  p <- 25
  N <- 50
  L <- matrix(rnorm(p * N), p, N)
  
  # Test symmetric pipeline
  aff_sym <- build_manifold_affinity(
    L, 
    k = 7,
    k_sigma = 10,  # Different from k
    normalization = "symmetric"
  )
  
  basis_sym <- compute_diffusion_basis(
    aff_sym,
    n_dims = 4,
    return_basis = TRUE,
    L_library_matrix = L
  )
  
  # Test row-stochastic pipeline
  aff_row <- build_manifold_affinity(
    L, 
    k = 7,
    normalization = "row_stochastic"
  )
  
  basis_row <- compute_diffusion_basis(
    aff_row,
    n_dims = 4,
    return_basis = TRUE,
    L_library_matrix = L
  )
  
  # Both should produce valid results
  expect_equal(dim(basis_sym$eigenvectors), c(N, 4))
  expect_equal(dim(basis_row$eigenvectors), c(N, 4))
  
  # Both should have 4 eigenvalues
  expect_equal(length(basis_sym$eigenvalues), 4)
  expect_equal(length(basis_row$eigenvalues), 4)
  
  # Eigenvalues should be different due to normalization
  expect_false(all(abs(basis_sym$eigenvalues - basis_row$eigenvalues) < 1e-6))
})

test_that("Memory efficiency: no dense matrices for large N", {
  skip_on_cran()  # Skip on CRAN due to memory
  
  set.seed(555)
  p <- 20
  N <- 1000  # Larger size
  
  # Generate in chunks to avoid memory issues
  L <- matrix(0, p, N)
  chunk_size <- 100
  for (i in seq(1, N, by = chunk_size)) {
    end_idx <- min(i + chunk_size - 1, N)
    L[, i:end_idx] <- matrix(rnorm(p * (end_idx - i + 1)), p)
  }
  
  # Should complete without memory explosion
  aff <- build_manifold_affinity(
    L, 
    k = 20,
    ann_method = "auto"
  )
  
  # Check it's sparse
  expect_s4_class(aff$W, "sparseMatrix")
  expect_s4_class(aff$T_norm, "sparseMatrix")
  
  # Sparsity should be high
  expect_gt(aff$diagnostics$sparsity, 0.95)
})

test_that("Numerical stability with identical points", {
  set.seed(666)
  p <- 10
  N <- 20
  
  # Create data with identical points to guarantee zero distances
  base_point <- matrix(rnorm(p), p, 1)
  L <- matrix(rep(base_point, N), p, N)  # All points identical
  
  # Should handle zero distances gracefully
  expect_warning(
    aff <- build_manifold_affinity(L, k = 5),
    "zero"
  )
  
  # Should still produce valid affinity
  expect_true(all(is.finite(aff$W@x)))
  expect_true(all(aff$W@x >= 0))
})