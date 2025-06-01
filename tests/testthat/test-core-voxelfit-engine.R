context("core voxel fit engine")

set.seed(1)

m <- 2
k <- 3
V <- 4
Gamma <- matrix(rnorm(m * k * V), m * k, V)

 test_that("extract_xi_beta_raw_svd_core returns correct dimensions", {
   res <- extract_xi_beta_raw_svd_core(Gamma, m, k)
   expect_equal(dim(res$Xi_raw_matrix), c(m, V))
   expect_equal(dim(res$Beta_raw_matrix), c(k, V))
 })

 test_that("apply_intrinsic_identifiability_core works", {
   Xi_raw <- matrix(rnorm(m * V), m, V)
   Beta_raw <- matrix(rnorm(k * V), k, V)
   B <- matrix(rnorm(5 * m), 5, m)
   h_ref <- rnorm(5)
   res <- apply_intrinsic_identifiability_core(Xi_raw, Beta_raw, B, h_ref)
   expect_equal(dim(res$Xi_ident_matrix), c(m, V))
   expect_equal(dim(res$Beta_ident_matrix), c(k, V))
 })

 test_that("make_voxel_graph_laplacian_core returns sparse Laplacian", {
   coords <- matrix(seq_len(9), ncol = 3)
   L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 2)
   expect_s4_class(L, "dgCMatrix")
   expect_equal(dim(L), c(nrow(coords), nrow(coords)))
 })

test_that("make_voxel_graph_laplacian_core produces correct Laplacian for 2-voxel chain", {
  coords <- matrix(c(0, 0, 0,
                     1, 0, 0), nrow = 2, byrow = TRUE)
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 1)
  Lm <- as.matrix(L)
  expect_equal(Lm, matrix(c(1, -1,
                             -1, 1), nrow = 2, byrow = TRUE))
})

test_that("apply_spatial_smoothing_core returns same dimension", {
  coords <- matrix(seq_len(9), ncol = 3)
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 2)
  Xi <- matrix(rnorm(m * nrow(coords)), m, nrow(coords))
  Xi_s <- apply_spatial_smoothing_core(Xi, L, 0.1)
  expect_equal(dim(Xi_s), dim(Xi))
})

test_that("prepare_lss_fixed_components_core returns matrices of correct size", {
  # Create a matrix where q_lss < n (3 columns, 10 rows)
  A <- cbind(1, matrix(rnorm(18), nrow = 10, ncol = 2))
  res <- prepare_lss_fixed_components_core(A, 1, 0.01)
  expect_equal(dim(res$P_lss_matrix), c(ncol(A), nrow(A)))
  expect_length(res$p_lss_vector, nrow(A))
})

test_that("reconstruct_hrf_shapes_core multiplies matrices", {
  B <- matrix(rnorm(10), 5, 2)
  Xi <- matrix(rnorm(2 * 3), 2, 3)
  H <- reconstruct_hrf_shapes_core(B, Xi)
  expect_equal(dim(H), c(5, 3))
})

test_that("run_lss_for_voxel_core returns vector of length T", {
  Y <- rnorm(5)
  X_list <- list(matrix(1:5, ncol = 1), matrix(5:1, ncol = 1))
  H <- rnorm(1)
  A <- cbind(1, matrix(rnorm(10), 5, 2))
  lss_fix <- prepare_lss_fixed_components_core(A, 1, 0.01)
  res <- run_lss_for_voxel_core(Y, X_list, H, A, lss_fix$P_lss_matrix, lss_fix$p_lss_vector)
  expect_length(res, length(X_list))
})

test_that("estimate_final_condition_betas_core returns matrix of correct dims", {
  Y <- matrix(rnorm(15), 5, 3)
  Xc <- list(matrix(1:5, ncol = 1), matrix(5:1, ncol = 1))
  H <- matrix(rnorm(1 * 3), 1, 3)
  res <- estimate_final_condition_betas_core(Y, Xc, H)
  expect_equal(dim(res), c(length(Xc), ncol(Y)))
})

