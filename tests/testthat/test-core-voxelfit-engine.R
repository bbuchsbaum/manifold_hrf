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

 test_that("apply_spatial_smoothing_core returns same dimension", {
   coords <- matrix(seq_len(9), ncol = 3)
   L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 2)
   Xi <- matrix(rnorm(m * nrow(coords)), m, nrow(coords))
   Xi_s <- apply_spatial_smoothing_core(Xi, L, 0.1)
   expect_equal(dim(Xi_s), dim(Xi))
 })

