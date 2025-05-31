# Tests for Core Manifold Construction Functions (Component 0)
# Tests for MHRF-CORE-MANIFOLD-01 and MHRF-CORE-MANIFOLD-02

test_that("calculate_manifold_affinity_core returns Markov matrix", {
  L <- matrix(seq(0, 4, length.out = 6), nrow = 2)
  S <- calculate_manifold_affinity_core(L, 1)
  expect_true(is.matrix(S))
  expect_equal(dim(S), c(ncol(L), ncol(L)))
})

test_that("get_manifold_basis_reconstructor_core returns reconstructor list", {
  L <- matrix(seq(0, 4, length.out = 6), nrow = 2)
  S <- calculate_manifold_affinity_core(L, 1)
  res <- get_manifold_basis_reconstructor_core(S, L, 2)
  expect_type(res, "list")
  expect_true(all(c("B_reconstructor_matrix", "Phi_coords_matrix") %in% names(res)))
  expect_equal(dim(res$B_reconstructor_matrix), c(nrow(L), 2))
})
