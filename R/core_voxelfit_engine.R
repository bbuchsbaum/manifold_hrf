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
  # Input validation
  if (!is.matrix(voxel_coords_matrix)) {
    stop("voxel_coords_matrix must be a matrix")
  }
  
  if (ncol(voxel_coords_matrix) != 3) {
    stop("voxel_coords_matrix must have exactly 3 columns (x, y, z coordinates)")
  }
  
  n_voxels <- nrow(voxel_coords_matrix)
  
  if (n_voxels < 2) {
    stop("voxel_coords_matrix must have at least 2 rows (voxels)")
  }
  
  if (!is.numeric(num_neighbors_Lsp) || length(num_neighbors_Lsp) != 1 || 
      num_neighbors_Lsp != round(num_neighbors_Lsp) || num_neighbors_Lsp < 1) {
    stop("num_neighbors_Lsp must be a positive integer")
  }
  
  # Handle edge case where we have fewer voxels than requested neighbors
  if (n_voxels <= num_neighbors_Lsp) {
    warning(sprintf("Requested %d neighbors but only %d other voxels available. Creating fully connected graph.",
                    num_neighbors_Lsp, n_voxels - 1))
    # Create fully connected graph
    W <- Matrix::Matrix(1, n_voxels, n_voxels) - Matrix::Diagonal(n_voxels)
  } else {
    nn <- RANN::nn2(voxel_coords_matrix, k = min(num_neighbors_Lsp + 1, n_voxels))
    i_idx <- rep(seq_len(n_voxels), each = ncol(nn$nn.idx) - 1)
    j_idx <- as.vector(nn$nn.idx[, -1])
    W <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1,
                              dims = c(n_voxels, n_voxels))
    W <- (W + Matrix::t(W)) / 2
  }
  
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
  # Input validation
  if (!is.matrix(Xi_ident_matrix)) {
    stop("Xi_ident_matrix must be a matrix")
  }
  
  if (!inherits(L_sp_sparse_matrix, c("Matrix", "sparseMatrix", "dgCMatrix"))) {
    stop("L_sp_sparse_matrix must be a sparse matrix (Matrix package)")
  }
  
  if (!is.numeric(lambda_spatial_smooth) || length(lambda_spatial_smooth) != 1 || 
      lambda_spatial_smooth < 0) {
    stop("lambda_spatial_smooth must be a non-negative scalar")
  }
  
  V <- ncol(Xi_ident_matrix)
  
  # Check dimensions match
  if (nrow(L_sp_sparse_matrix) != V || ncol(L_sp_sparse_matrix) != V) {
    stop(sprintf("L_sp_sparse_matrix must be %d x %d to match Xi_ident_matrix with %d voxels",
                 V, V, V))
  }
  
  A <- Matrix::Diagonal(V) + lambda_spatial_smooth * L_sp_sparse_matrix
  Xi_smoothed <- Xi_ident_matrix
  for (j in seq_len(nrow(Xi_ident_matrix))) {
    Xi_smoothed[j, ] <- as.vector(Matrix::solve(A, Xi_ident_matrix[j, ]))
  }
  Xi_smoothed
}

#' Precompute fixed components for Woodbury LSS kernel
#'
#' @param A_lss_fixed_matrix n x q_lss matrix of nuisance regressors
#' @param intercept_col_index_in_Alss integer index of intercept column or NULL
#' @param lambda_ridge_Alss ridge penalty for inversion stability
#' @return list with P_lss_matrix (q_lss x n) and p_lss_vector (n length)
#' @export
prepare_lss_fixed_components_core <- function(A_lss_fixed_matrix,
                                              intercept_col_index_in_Alss = NULL,
                                              lambda_ridge_Alss = 0) {
  AtA <- crossprod(A_lss_fixed_matrix)
  if (lambda_ridge_Alss > 0) {
    AtA <- AtA + diag(lambda_ridge_Alss, ncol(A_lss_fixed_matrix))
  }
  jitter <- 1e-6 * median(diag(AtA))
  AtA <- AtA + diag(jitter, ncol(A_lss_fixed_matrix))
  P_lss <- solve(AtA, t(A_lss_fixed_matrix))
  if (!is.null(intercept_col_index_in_Alss)) {
    ginvA <- MASS::ginv(A_lss_fixed_matrix)
    p_vec <- ginvA[intercept_col_index_in_Alss, ]
  } else {
    p_vec <- rep(0, nrow(A_lss_fixed_matrix))
  }
  list(P_lss_matrix = P_lss, p_lss_vector = p_vec)
}

#' Reconstruct voxel-specific HRF shapes from manifold coordinates
#'
#' @param B_reconstructor_matrix p x m manifold reconstructor
#' @param Xi_smoothed_matrix m x V matrix of smoothed manifold coordinates
#' @return p x V matrix of HRF shapes
#' @export
reconstruct_hrf_shapes_core <- function(B_reconstructor_matrix,
                                        Xi_smoothed_matrix) {
  B_reconstructor_matrix %*% Xi_smoothed_matrix
}

#' Woodbury LSS kernel for a single voxel
#'
#' @param Y_proj_voxel_vector n-length projected BOLD vector
#' @param X_trial_onset_list_of_matrices list of T n x p trial onset matrices
#' @param H_shape_voxel_vector p-length HRF shape vector for voxel
#' @param A_lss_fixed_matrix n x q_lss fixed nuisance regressors
#' @param P_lss_matrix q_lss x n matrix from prepare_lss_fixed_components_core
#' @param p_lss_vector n-length vector from prepare_lss_fixed_components_core
#' @return T-length vector of trial betas for this voxel
#' @export
run_lss_for_voxel_core <- function(Y_proj_voxel_vector,
                                   X_trial_onset_list_of_matrices,
                                   H_shape_voxel_vector,
                                   A_lss_fixed_matrix,
                                   P_lss_matrix,
                                   p_lss_vector) {
  T_trials <- length(X_trial_onset_list_of_matrices)
  n <- length(Y_proj_voxel_vector)
  C_v <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shape_voxel_vector
  }
  U_v <- P_lss_matrix %*% C_v
  V_reg <- C_v - A_lss_fixed_matrix %*% U_v
  pc_row <- drop(p_lss_vector %*% C_v)
  cv_row <- colSums(V_reg * V_reg)
  alpha_row <- (1 - pc_row) / pmax(cv_row, .Machine$double.eps)
  S_eff <- sweep(V_reg, 2, alpha_row, "*")
  S_eff <- sweep(S_eff, 1, p_lss_vector, "+")
  drop(crossprod(S_eff, Y_proj_voxel_vector))
}

#' Run LSS voxel loop over all voxels
#'
#' @param Y_proj_matrix n x V projected BOLD matrix
#' @param X_trial_onset_list_of_matrices list of T n x p trial onset matrices
#' @param H_shapes_allvox_matrix p x V matrix of voxel HRF shapes
#' @param A_lss_fixed_matrix n x q_lss fixed nuisance regressors
#' @param P_lss_matrix q_lss x n matrix precomputed
#' @param p_lss_vector n-length vector precomputed
#' @param ram_heuristic_GB_for_Rt numeric, ignored in simplified implementation
#' @return T x V matrix of trial betas
#' @export
run_lss_voxel_loop_core <- function(Y_proj_matrix,
                                    X_trial_onset_list_of_matrices,
                                    H_shapes_allvox_matrix,
                                    A_lss_fixed_matrix,
                                    P_lss_matrix,
                                    p_lss_vector,
                                    ram_heuristic_GB_for_Rt = Inf) {
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  Beta_trial_allvox <- matrix(0, T_trials, V)
  for (v in seq_len(V)) {
    Beta_trial_allvox[, v] <- run_lss_for_voxel_core(
      Y_proj_matrix[, v],
      X_trial_onset_list_of_matrices,
      H_shapes_allvox_matrix[, v],
      A_lss_fixed_matrix,
      P_lss_matrix,
      p_lss_vector
    )
  }
  Beta_trial_allvox
}

#' Estimate final condition-level betas using smoothed HRFs
#'
#' @param Y_proj_matrix n x V projected BOLD matrix
#' @param X_condition_list_proj_matrices list of k n x p design matrices
#' @param H_shapes_allvox_matrix p x V HRF shapes for all voxels
#' @param lambda_beta_final ridge penalty
#' @param control_alt_list list with max_iter (ignored here)
#' @return k x V matrix of final betas
#' @export
estimate_final_condition_betas_core <- function(Y_proj_matrix,
                                                X_condition_list_proj_matrices,
                                                H_shapes_allvox_matrix,
                                                lambda_beta_final = 0,
                                                control_alt_list = list(max_iter = 1)) {
  # Input validation
  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }
  if (!is.list(X_condition_list_proj_matrices)) {
    stop("X_condition_list_proj_matrices must be a list")
  }
  if (length(X_condition_list_proj_matrices) == 0) {
    stop("X_condition_list_proj_matrices must contain at least one condition")
  }
  if (!is.matrix(H_shapes_allvox_matrix)) {
    stop("H_shapes_allvox_matrix must be a matrix")
  }
  
  k <- length(X_condition_list_proj_matrices)
  V <- ncol(Y_proj_matrix)
  n <- nrow(Y_proj_matrix)
  
  # More validation
  if (ncol(H_shapes_allvox_matrix) != V) {
    stop("H_shapes_allvox_matrix must have V columns matching Y_proj_matrix")
  }
  if (!is.numeric(lambda_beta_final) || length(lambda_beta_final) != 1 || lambda_beta_final < 0) {
    stop("lambda_beta_final must be a non-negative scalar")
  }
  if (!is.list(control_alt_list)) {
    stop("control_alt_list must be a list")
  }
  if (!is.null(control_alt_list$max_iter)) {
    if (!is.numeric(control_alt_list$max_iter) || control_alt_list$max_iter < 1) {
      stop("max_iter must be a positive integer")
    }
  }
  if (!is.null(control_alt_list$rel_change_tol)) {
    if (!is.numeric(control_alt_list$rel_change_tol) || control_alt_list$rel_change_tol < 0) {
      stop("rel_change_tol must be a non-negative scalar")
    }
  }
  
  Beta_final <- matrix(0, k, V)
  for (v in seq_len(V)) {
    X_v <- matrix(0, n, k)
    for (c in seq_len(k)) {
      X_v[, c] <- X_condition_list_proj_matrices[[c]] %*% H_shapes_allvox_matrix[, v]
    }
    XtX <- crossprod(X_v) + diag(lambda_beta_final, k)
    XtY <- crossprod(X_v, Y_proj_matrix[, v])
    
    # Check for near-singularity
    if (rcond(XtX) < .Machine$double.eps) {
      # If singular, return zeros
      Beta_final[, v] <- rep(0, k)
    } else {
      Beta_final[, v] <- solve(XtX, XtY)
    }
  }
  Beta_final
}

