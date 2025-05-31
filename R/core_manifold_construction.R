# Core Manifold Construction Functions (Component 0)
# Implementation of MHRF-CORE-MANIFOLD-01 and MHRF-CORE-MANIFOLD-02

#' Calculate Manifold Affinity and Markov Matrix (Core)
#'
#' @param L_library_matrix A p x N matrix of HRF shapes
#' @param k_local_nn_for_sigma Integer, k-nearest neighbors for self-tuning bandwidth
#' @param use_sparse_W_params List with parameters for sparse W matrix
#' @return S_markov_matrix (N x N matrix or sparse matrix)
#' @export
calculate_manifold_affinity_core <- function(L_library_matrix,
                                            k_local_nn_for_sigma,
                                            use_sparse_W_params = list()) {
  if (!is.matrix(L_library_matrix)) {
    stop("L_library_matrix must be a matrix")
  }
  if (k_local_nn_for_sigma <= 0) {
    stop("k_local_nn_for_sigma must be positive")
  }

  N <- ncol(L_library_matrix)
  # pairwise Euclidean distances between HRF shapes
  dist_mat <- as.matrix(dist(t(L_library_matrix), method = "euclidean"))

  # self tuning bandwidth sigma_i based on k-th nearest neighbour
  dist_sorted <- apply(dist_mat, 1, sort)
  sigma <- dist_sorted[k_local_nn_for_sigma + 1, ]

  W <- matrix(0, N, N)
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      W[i, j] <- exp(- (dist_mat[i, j]^2) / (sigma[i] * sigma[j]))
    }
  }

  D_inv <- diag(1 / rowSums(W))
  S <- D_inv %*% W
  S
}

#' Get Manifold Basis Reconstructor (Core)
#'
#' @param S_markov_matrix N x N Markov matrix
#' @param L_library_matrix p x N matrix of HRF shapes
#' @param m_manifold_dim_target Target manifold dimensionality
#' @param m_manifold_dim_min_variance Minimum variance explained threshold
#' @return List with B_reconstructor_matrix, Phi_coords_matrix, eigenvalues, etc.
#' @export
get_manifold_basis_reconstructor_core <- function(S_markov_matrix,
                                                 L_library_matrix,
                                                 m_manifold_dim_target,
                                                 m_manifold_dim_min_variance = 0.95) {
  if (!is.matrix(S_markov_matrix) || !is.matrix(L_library_matrix)) {
    stop("Inputs must be matrices")
  }

  eig <- eigen(S_markov_matrix)
  evals <- Re(eig$values)
  evecs <- Re(eig$vectors)
  ord <- order(evals, decreasing = TRUE)
  evals <- evals[ord]
  evecs <- evecs[, ord, drop = FALSE]

  cum_var <- cumsum(evals[-1]) / sum(evals[-1])
  m_auto <- which(cum_var >= m_manifold_dim_min_variance)[1]
  if (is.na(m_auto)) {
    m_auto <- length(cum_var)
  }
  m_final <- min(m_manifold_dim_target, m_auto)
  Phi <- evecs[, 2:(m_final + 1), drop = FALSE]

  B_rec <- L_library_matrix %*% Phi %*%
    solve(crossprod(Phi) + 1e-8 * diag(m_final))

  list(B_reconstructor_matrix = B_rec,
       Phi_coords_matrix = Phi,
       eigenvalues_S_vector = evals,
       m_final_dim = m_final,
       m_auto_selected_dim = m_auto)
}
