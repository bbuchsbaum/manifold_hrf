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

