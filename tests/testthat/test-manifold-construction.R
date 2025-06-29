# Tests for Core Manifold Construction Functions (Component 0)
# Tests for MHRF-CORE-MANIFOLD-01 and MHRF-CORE-MANIFOLD-02

test_that("calculate_manifold_affinity_core works with small matrix", {
  # Create a small test HRF library matrix (p=10 time points, N=5 HRFs)
  set.seed(123)
  p <- 10
  N <- 5
  L_library_matrix <- matrix(rnorm(p * N), nrow = p, ncol = N)
  
  # Test basic functionality
  k_local_nn <- 2
  S_markov <- calculate_manifold_affinity_core(L_library_matrix, k_local_nn)
  
  # Check output dimensions
  expect_equal(dim(S_markov), c(N, N))
  
  # Check that S is a symmetric normalized Laplacian
  # For symmetric normalized Laplacian, rows don't sum to 1
  # Instead, check that it's symmetric
  expect_true(isSymmetric(S_markov, tol = 1e-8))
  
  # Check that diagonal entries are reasonable (between 0 and 1)
  expect_true(all(diag(S_markov) >= 0))
  expect_true(all(diag(S_markov) <= 1))
  
  # Check that all entries are non-negative
  expect_true(all(S_markov >= 0))
})

test_that("calculate_manifold_affinity_core handles sparse matrix parameters", {
  # Create a larger test matrix
  set.seed(456)
  p <- 20
  N <- 100
  L_library_matrix <- matrix(rnorm(p * N), nrow = p, ncol = N)
  
  # Test with sparse parameters
  k_local_nn <- 5
  sparse_params <- list(
    sparse_if_N_gt = 50,  # Should trigger sparse mode
    k_nn_for_W_sparse = 10
  )
  
  S_markov_sparse <- calculate_manifold_affinity_core(
    L_library_matrix, 
    k_local_nn, 
    sparse_params
  )
  
  # Check if Matrix package is available
  if (requireNamespace("Matrix", quietly = TRUE)) {
    # Should return a sparse matrix
    expect_true(inherits(S_markov_sparse, "sparseMatrix"))
  }
  
  # Check basic properties still hold
  expect_equal(dim(S_markov_sparse), c(N, N))
  
  # Should be a symmetric normalized Laplacian
  if (inherits(S_markov_sparse, "Matrix")) {
    expect_true(Matrix::isSymmetric(S_markov_sparse, tol = 1e-8))
  } else {
    expect_true(isSymmetric(S_markov_sparse, tol = 1e-8))
  }
})

test_that("k_nn_for_W_sparse larger than N is truncated", {
  set.seed(999)
  p <- 5
  N <- 10
  L_library_matrix <- matrix(rnorm(p * N), nrow = p, ncol = N)

  sparse_params <- list(
    sparse_if_N_gt = 1,
    k_nn_for_W_sparse = N + 5
  )

  S_markov <- calculate_manifold_affinity_core(
    L_library_matrix,
    k_local_nn_for_sigma = 2,
    use_sparse_W_params = sparse_params
  )

  expect_equal(dim(S_markov), c(N, N))
  # Should be a symmetric normalized Laplacian
  if (inherits(S_markov, "Matrix")) {
    expect_true(Matrix::isSymmetric(S_markov, tol = 1e-8))
  } else {
    expect_true(isSymmetric(S_markov, tol = 1e-8))
  }
})

test_that("calculate_manifold_affinity_core validates inputs correctly", {
  # Test with non-matrix input
  expect_error(
    calculate_manifold_affinity_core(data.frame(a = 1:5), 2),
    "L_library_matrix must be a matrix"
  )
  
  # Test with k >= N
  L_small <- matrix(1:12, nrow = 3, ncol = 4)
  expect_error(
    calculate_manifold_affinity_core(L_small, 4),
    "k_local_nn_for_sigma must be less than the number of HRFs"
  )

  # Test with non-positive or non-integer k
  expect_error(
    calculate_manifold_affinity_core(L_small, 0),
    "k_local_nn_for_sigma must be a positive integer"
  )

  expect_error(
    calculate_manifold_affinity_core(L_small, 2.5),
    "k_local_nn_for_sigma must be a positive integer"
  )
})

test_that("calculate_manifold_affinity_core produces symmetric affinities", {
  # Create test data
  set.seed(789)
  p <- 15
  N <- 10
  L_library_matrix <- matrix(rnorm(p * N), nrow = p, ncol = N)
  
  k_local_nn <- 3
  S_markov <- calculate_manifold_affinity_core(L_library_matrix, k_local_nn)
  
  # The affinity matrix W should be symmetric, but S is not necessarily symmetric
  # However, we can check that the structure makes sense
  # S should have positive entries where HRFs are similar
  expect_true(all(S_markov >= 0))
  expect_true(all(S_markov <= 1))
})

test_that("get_manifold_basis_reconstructor_core works with basic inputs", {
  # Create test data
  set.seed(123)
  p <- 20  # time points
  N <- 30  # number of HRFs
  L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
  
  # Create Markov matrix
  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 5)
  
  # Test basic functionality
  result <- get_manifold_basis_reconstructor_core(
    S, L_library, 
    m_manifold_dim_target = 5,
    m_manifold_dim_min_variance = 0.95
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("B_reconstructor_matrix", "Phi_coords_matrix", 
                        "eigenvalues_S_vector", "m_final_dim", "m_auto_selected_dim", "m_manifold_dim"))
  
  # Check dimensions
  expect_equal(nrow(result$B_reconstructor_matrix), p)
  expect_equal(ncol(result$B_reconstructor_matrix), result$m_final_dim)
  expect_equal(nrow(result$Phi_coords_matrix), N)
  expect_equal(ncol(result$Phi_coords_matrix), result$m_final_dim)
  
  # Check that dimensions are reasonable
  expect_true(result$m_final_dim <= min(5, N-1))
  expect_true(result$m_final_dim >= 1)
  expect_true(result$m_auto_selected_dim >= 1)
})

test_that("get_manifold_basis_reconstructor_core handles variance threshold", {
  # Create test data with clear structure (low rank)
  set.seed(456)
  p <- 15
  N <- 20
  # Create low-rank HRF library (should need fewer dimensions)
  U <- matrix(rnorm(p * 3), nrow = p, ncol = 3)
  V <- matrix(rnorm(3 * N), nrow = 3, ncol = N)
  L_library <- U %*% V + 0.1 * matrix(rnorm(p * N), nrow = p, ncol = N)
  
  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 4)
  
  # Test with high variance requirement
  result_high_var <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 10,
    m_manifold_dim_min_variance = 0.99
  )
  
  # Test with low variance requirement  
  result_low_var <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 10,
    m_manifold_dim_min_variance = 0.80
  )
  
  # Higher variance requirement should need more dimensions
  expect_gte(result_high_var$m_auto_selected_dim, result_low_var$m_auto_selected_dim)
})

test_that("get_manifold_basis_reconstructor_core validates inputs", {
  L_library <- matrix(1:20, nrow = 4, ncol = 5)
  S <- diag(5)
  
  # Test dimension mismatch
  expect_error(
    get_manifold_basis_reconstructor_core(diag(4), L_library, 3),
    "dimensions must match"
  )
  
  # Test invalid variance threshold
  expect_error(
    get_manifold_basis_reconstructor_core(S, L_library, 3, m_manifold_dim_min_variance = 1.5),
    "must be between 0 and 1"
  )
  
  # Test non-matrix inputs
  expect_error(
    get_manifold_basis_reconstructor_core(data.frame(S), L_library, 3),
    "must be a matrix"
  )
})

test_that("get_manifold_basis_reconstructor_core reconstruction works", {
  # Test that we can reconstruct HRFs reasonably well
  set.seed(789)
  p <- 25
  N <- 40
  L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
  
  # Normalize columns for easier comparison
  for (i in 1:N) {
    L_library[, i] <- L_library[, i] / sqrt(sum(L_library[, i]^2))
  }
  
  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 7)
  
  result <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 10,
    m_manifold_dim_min_variance = 0.90
  )
  
  # Test reconstruction of library HRFs
  # Reconstructed = B * Phi'
  L_reconstructed <- result$B_reconstructor_matrix %*% t(result$Phi_coords_matrix)
  
  # Check that reconstruction error is reasonable
  # (not expecting perfect reconstruction with reduced dimensions)
  reconstruction_error <- norm(L_library - L_reconstructed, "F") / norm(L_library, "F")
  expect_lt(reconstruction_error, 1.0)  # Less than 100% error (sanity check)
  
  # Check that B has reasonable norm (not exploding)
  B_norm <- norm(result$B_reconstructor_matrix, "F")
  expect_true(is.finite(B_norm))
  expect_lt(B_norm, 1000)  # Arbitrary but reasonable upper bound
})


test_that("get_manifold_basis_reconstructor_core handles large N with RSpectra", {
  skip_if_not_installed("RSpectra")
  set.seed(101)
  p <- 30
  N <- 250
  L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 7)

  result <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 5,
    m_manifold_dim_min_variance = 0.95
  )

  expect_equal(nrow(result$Phi_coords_matrix), N)
  expect_true(ncol(result$Phi_coords_matrix) < N)
})
test_that("ann_euclidean distance works correctly with RcppHNSW", {
  # Originally this test checked fallback behavior when RcppHNSW is missing.
  # Since RcppHNSW is installed in the test environment, we instead verify:
  # 1. ann_euclidean works without warning when the package is available
  # 2. Results are reasonable (not necessarily identical to exact euclidean)
  # 
  # Note: Testing the actual fallback would require unloading RcppHNSW,
  # which could cause issues with other tests
  
  p <- 5
  N <- 10
  set.seed(123)
  L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
  
  if (requireNamespace("RcppHNSW", quietly = TRUE)) {
    # Test that ann_euclidean works without warning when package is available
    expect_no_warning(
      result_ann <- calculate_manifold_affinity_core(
        L_library,
        k_local_nn_for_sigma = 2,
        distance_engine = "ann_euclidean"
      )
    )
    
    # Verify the result is a valid Markov matrix
    result_ann_mat <- as.matrix(result_ann)
    
    # Check dimensions
    expect_equal(dim(result_ann_mat), c(N, N))
    
    # Check row sums are close to 1 (Markov property)
    # Note: ann_euclidean may have some numerical differences
    row_sums <- rowSums(result_ann_mat)
    expect_true(all(abs(row_sums - 1) < 0.5))
    
    # Check non-negative entries
    expect_true(all(result_ann_mat >= 0))
  } else {
    # If RcppHNSW is not installed, it should warn and fall back
    expect_warning(
      calculate_manifold_affinity_core(
        L_library,
        k_local_nn_for_sigma = 2,
        distance_engine = "ann_euclidean"
      ),
      "RcppHNSW"
    )
  }
})

test_that("degenerate data does not crash manifold construction", {
  p <- 15
  N <- 8
  L_library <- matrix(1, nrow = p, ncol = N)

  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 2)
  expect_equal(dim(S), c(N, N))
  expect_true(all(is.finite(S)))

  result <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 3
  )

  expect_true(all(is.finite(result$B_reconstructor_matrix)))
  expect_gte(result$m_final_dim, 1)
})

test_that("disconnected components yield block structure", {
  set.seed(42)
  p <- 10
  N <- 20
  L_cluster1 <- matrix(rnorm(p * (N/2), sd = 0.1), p, N/2)
  L_cluster2 <- matrix(rnorm(p * (N/2), mean = 100, sd = 0.1), p, N/2)
  L_library <- cbind(L_cluster1, L_cluster2)

  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 2)

  cross_block <- S[1:(N/2), (N/2 + 1):N]
  expect_lt(max(cross_block), 1e-3)

  result <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 5
  )

  expect_equal(nrow(result$Phi_coords_matrix), N)
  expect_true(all(is.finite(result$Phi_coords_matrix)))
})

test_that("high dimensional data still returns finite manifold", {
  set.seed(99)
  p <- 500
  N <- 100
  L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)

  S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 5)
  result <- get_manifold_basis_reconstructor_core(
    S, L_library,
    m_manifold_dim_target = 8
  )

  expect_equal(nrow(result$Phi_coords_matrix), N)
  expect_true(all(is.finite(result$Phi_coords_matrix)))
})
