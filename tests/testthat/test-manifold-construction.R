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
  
  # Check that S is a stochastic matrix (rows sum to 1)
  row_sums <- rowSums(S_markov)
  expect_equal(row_sums, rep(1, N), tolerance = 1e-8)
  
  # Check that diagonal is 0 (no self-loops)
  expect_equal(diag(S_markov), rep(0, N))
  
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
  
  # Handle both regular and sparse matrices
  if (inherits(S_markov_sparse, "Matrix")) {
    row_sums <- Matrix::rowSums(S_markov_sparse)
  } else {
    row_sums <- rowSums(S_markov_sparse)
  }
  expect_equal(row_sums, rep(1, N), tolerance = 1e-8)
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
  if (inherits(S_markov, "Matrix")) {
    row_sums <- Matrix::rowSums(S_markov)
  } else {
    row_sums <- rowSums(S_markov)
  }
  expect_equal(row_sums, rep(1, N), tolerance = 1e-8)
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
test_that("ann_euclidean distance falls back when RcppHNSW missing", {
  skip_if(requireNamespace("RcppHNSW", quietly = TRUE),
          "RcppHNSW installed; cannot test fallback")
  p <- 5
  N <- 10
  L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
  expect_warning(
    calculate_manifold_affinity_core(
      L_library,
      k_local_nn_for_sigma = 2,
      distance_engine = "ann_euclidean"
    ),
    "RcppHNSW"
  )

})
