# Core Voxel-wise Fit Engine Functions (Component 1)

#' Project out confound regressors from data and design matrices
#'
#' @param Y_data_matrix n x V numeric matrix of BOLD data
#' @param X_list_of_matrices list of design matrices (each n x p)
#' @param Z_confounds_matrix optional n x q matrix of confounds
#' @return list with Y_proj_matrix and X_list_proj_matrices
#' @export
project_out_confounds_core <- function(Y_data_matrix,
                                       X_list_of_matrices,
                                       Z_confounds_matrix = NULL) {
  if (is.null(Z_confounds_matrix)) {
    return(list(Y_proj_matrix = Y_data_matrix,
                X_list_proj_matrices = X_list_of_matrices))
  }
  Qz <- qr.Q(qr(Z_confounds_matrix))
  Y_proj <- Y_data_matrix - Qz %*% (t(Qz) %*% Y_data_matrix)
  X_proj <- lapply(X_list_of_matrices, function(X) {
    X - Qz %*% (t(Qz) %*% X)
  })
  list(Y_proj_matrix = Y_proj,
       X_list_proj_matrices = X_proj)
}

#' Transform design matrices to manifold basis
#'
#' @param X_condition_list_proj_matrices list of projected design matrices (n x p)
#' @param B_reconstructor_matrix p x m manifold reconstructor matrix
#' @return list of design matrices in manifold basis
#' @export
transform_designs_to_manifold_basis_core <- function(X_condition_list_proj_matrices,
                                                     B_reconstructor_matrix) {
  lapply(X_condition_list_proj_matrices, function(X) {
    X %*% B_reconstructor_matrix
  })
}

#' Solve GLM for gamma coefficients
#'
#' @param Z_list_of_matrices list of n x m design matrices in manifold basis
#' @param Y_proj_matrix n x V projected BOLD matrix
#' @param lambda_gamma ridge penalty
#' @param orthogonal_approx_flag logical indicating orthogonal design approximation
#' @return (k*m) x V matrix of gamma coefficients
#' @export
solve_glm_for_gamma_core <- function(Z_list_of_matrices,
                                     Y_proj_matrix,
                                     lambda_gamma = 0,
                                     orthogonal_approx_flag = FALSE) {
  k <- length(Z_list_of_matrices)
  n <- nrow(Y_proj_matrix)
  m <- ncol(Z_list_of_matrices[[1]])
  Xt <- do.call(cbind, Z_list_of_matrices)
  XtX <- crossprod(Xt)
  if (lambda_gamma > 0) {
    XtX <- XtX + diag(lambda_gamma, k * m)
  }
  XtY <- crossprod(Xt, Y_proj_matrix)
  beta <- solve(XtX, XtY)
  beta
}


#' Extract raw manifold coordinates and condition amplitudes via SVD
#'
#' @param Gamma_coeffs_matrix (k*m) x V matrix of gamma coefficients
#' @param m_manifold_dim Integer manifold dimension m
#' @param k_conditions Integer number of conditions k
#' @return list with Xi_raw_matrix (m x V) and Beta_raw_matrix (k x V)
#' @export
extract_xi_beta_raw_svd_core <- function(Gamma_coeffs_matrix,
                                         m_manifold_dim,
                                         k_conditions) {
  if (!is.matrix(Gamma_coeffs_matrix)) {
    stop("Gamma_coeffs_matrix must be a matrix")
  }
  if (nrow(Gamma_coeffs_matrix) != m_manifold_dim * k_conditions) {
    stop("nrow(Gamma_coeffs_matrix) must equal m*k")
  }
  V <- ncol(Gamma_coeffs_matrix)
  Xi_raw <- matrix(0, m_manifold_dim, V)
  Beta_raw <- matrix(0, k_conditions, V)
  for (v in seq_len(V)) {
    Gv <- matrix(Gamma_coeffs_matrix[, v], nrow = m_manifold_dim, ncol = k_conditions)
    sv <- svd(Gv)
    if (length(sv$d) == 0 || sv$d[1] < .Machine$double.eps) {
      next
    }
    Xi_raw[, v] <- sv$u[, 1] * sqrt(sv$d[1])
    Beta_raw[, v] <- sv$v[, 1] * sqrt(sv$d[1])
  }
  list(Xi_raw_matrix = Xi_raw, Beta_raw_matrix = Beta_raw)
}

#' Apply intrinsic identifiability constraints
#'
#' @param Xi_raw_matrix m x V matrix of raw manifold coordinates
#' @param Beta_raw_matrix k x V matrix of raw condition amplitudes
#' @param B_reconstructor_matrix p x m manifold reconstructor
#' @param h_ref_shape_vector p-length canonical HRF shape
#' @param ident_scale_method one of "l2_norm", "max_abs_val", "none"
#' @param ident_sign_method one of "first_component", "canonical_correlation"
#' @return list with Xi_ident_matrix and Beta_ident_matrix
#' @export
apply_intrinsic_identifiability_core <- function(Xi_raw_matrix,
                                                 Beta_raw_matrix,
                                                 B_reconstructor_matrix,
                                                 h_ref_shape_vector,
                                                 ident_scale_method = c("l2_norm", "max_abs_val", "none"),
                                                 ident_sign_method = c("first_component", "canonical_correlation")) {
  ident_scale_method <- match.arg(ident_scale_method)
  ident_sign_method <- match.arg(ident_sign_method)

  xi_ref_coord <- MASS::ginv(B_reconstructor_matrix) %*% h_ref_shape_vector

  V <- ncol(Xi_raw_matrix)
  Xi_ident <- matrix(0, nrow(Xi_raw_matrix), V)
  Beta_ident <- matrix(0, nrow(Beta_raw_matrix), V)

  for (v in seq_len(V)) {
    xi_v <- Xi_raw_matrix[, v]
    beta_v <- Beta_raw_matrix[, v]

    if (ident_sign_method == "canonical_correlation") {
      sgn <- sign(sum(xi_v * xi_ref_coord))
      if (sgn == 0) sgn <- 1
    } else { # first_component
      sgn <- sign(xi_v[1])
      if (sgn == 0) sgn <- 1
    }

    xi_v <- xi_v * sgn
    beta_v <- beta_v * sgn

    hrf_v <- B_reconstructor_matrix %*% xi_v
    scale_val <- 1
    if (ident_scale_method == "l2_norm") {
      scale_val <- 1 / pmax(sqrt(sum(hrf_v^2)), .Machine$double.eps)
    } else if (ident_scale_method == "max_abs_val") {
      scale_val <- 1 / pmax(max(abs(hrf_v)), .Machine$double.eps)
    }
    Xi_ident[, v] <- xi_v * scale_val
    Beta_ident[, v] <- beta_v / scale_val
  }

  list(Xi_ident_matrix = Xi_ident, Beta_ident_matrix = Beta_ident)
}

#' Construct voxel graph Laplacian
#'
#' @param voxel_coords_matrix V x 3 matrix of voxel coordinates
#' @param num_neighbors_Lsp number of nearest neighbours
#' @return sparse V x V graph Laplacian matrix
#' @export
make_voxel_graph_laplacian_core <- function(voxel_coords_matrix, num_neighbors_Lsp = 6) {
  nn <- RANN::nn2(voxel_coords_matrix, k = num_neighbors_Lsp + 1)
  i_idx <- rep(seq_len(nrow(voxel_coords_matrix)), each = num_neighbors_Lsp)
  j_idx <- as.vector(nn$nn.idx[, -1])
  W <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1,
                            dims = c(nrow(voxel_coords_matrix), nrow(voxel_coords_matrix)))
  W <- (W + t(W)) / 2
  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  L <- D - W
  L
}

#' Apply spatial smoothing to manifold coordinates
#'
#' @param Xi_ident_matrix m x V matrix of manifold coordinates
#' @param L_sp_sparse_matrix V x V Laplacian matrix
#' @param lambda_spatial_smooth smoothing strength
#' @return Xi_smoothed_matrix m x V matrix
#' @export
apply_spatial_smoothing_core <- function(Xi_ident_matrix,
                                         L_sp_sparse_matrix,
                                         lambda_spatial_smooth) {
  V <- ncol(Xi_ident_matrix)
  A <- Matrix::Diagonal(V) + lambda_spatial_smooth * L_sp_sparse_matrix
  Xi_smoothed <- Xi_ident_matrix
  for (j in seq_len(nrow(Xi_ident_matrix))) {
    Xi_smoothed[j, ] <- as.vector(Matrix::solve(A, Xi_ident_matrix[j, ]))
  }
  Xi_smoothed
}

