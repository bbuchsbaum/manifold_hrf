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

  if (!is.matrix(Y_data_matrix)) {
    stop("Y_data_matrix must be a matrix")
  }

  if (!is.list(X_list_of_matrices)) {
    stop("X_list_of_matrices must be a list")
  }

  n <- nrow(Y_data_matrix)

  for (i in seq_along(X_list_of_matrices)) {
    if (!is.matrix(X_list_of_matrices[[i]]) || nrow(X_list_of_matrices[[i]]) != n) {
      stop("X_list_of_matrices[[", i, "]] must be a matrix with ", n, " rows")
    }
  }

  if (is.null(Z_confounds_matrix)) {
    return(list(Y_proj_matrix = Y_data_matrix,
                X_list_proj_matrices = X_list_of_matrices))
  }

  if (!is.matrix(Z_confounds_matrix)) {
    stop("Z_confounds_matrix must be a matrix or NULL")
  }

  if (anyNA(Z_confounds_matrix)) {
    stop("Z_confounds_matrix must not contain NA values")
  }

  if (nrow(Z_confounds_matrix) != n) {
    stop("Z_confounds_matrix must have the same number of rows as Y_data_matrix")
  }

  if (ncol(Z_confounds_matrix) >= n) {
    stop("Z_confounds_matrix has too many columns (must be less than number of timepoints)")
  }

  qr_Z <- qr(Z_confounds_matrix, LAPACK = TRUE)
  if (qr_Z$rank < ncol(Z_confounds_matrix)) {
    warning("Z_confounds_matrix is rank deficient; using independent columns only")
  }
  Qz <- qr.Q(qr_Z)[, seq_len(qr_Z$rank), drop = FALSE]
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

  if (!is.list(X_condition_list_proj_matrices)) {
    stop("X_condition_list_proj_matrices must be a list")
  }

  if (length(X_condition_list_proj_matrices) == 0) {
    stop("X_condition_list_proj_matrices cannot be empty")
  }

  if (!is.matrix(B_reconstructor_matrix)) {
    stop("B_reconstructor_matrix must be a matrix")
  }

  p <- nrow(B_reconstructor_matrix)
  m <- ncol(B_reconstructor_matrix)

  n <- nrow(X_condition_list_proj_matrices[[1]])

  for (i in seq_along(X_condition_list_proj_matrices)) {
    X <- X_condition_list_proj_matrices[[i]]
    if (!is.matrix(X)) {
      stop(sprintf("X_condition_list_proj_matrices[[%d]] must be a matrix", i))
    }
    if (nrow(X) != n) {
      stop(sprintf("All X matrices must have %d rows", n))
    }
    if (ncol(X) != p) {
      stop(sprintf("X_condition_list_proj_matrices[[%d]] has %d columns but B_reconstructor_matrix has %d rows",
                   i, ncol(X), p))
    }
  }

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

#' Robust SVD Extraction with Conditioning
#'
#' Wrapper for extracting Xi and Beta using a numerically stable SVD.
#' This version includes automatic regularization and fallback strategies.
#'
#' @param Gamma_coeffs_matrix (k*m) x V matrix of gamma coefficients
#' @param m_manifold_dim Integer manifold dimension m
#' @param k_conditions Integer number of conditions k
#' @param regularization_factor Multiplier for diagonal regularization
#' @param max_condition_number Threshold for conditioning warning
#' @param use_randomized_svd Logical, use randomized SVD if available
#' @param logger Optional logger object
#' @return list with Xi_raw_matrix, Beta_raw_matrix, and quality metrics
#' @export
extract_xi_beta_raw_svd_robust <- function(Gamma_coeffs_matrix,
                                           m_manifold_dim,
                                           k_conditions,
                                           regularization_factor = 10,
                                           max_condition_number = 1e8,
                                           use_randomized_svd = FALSE,
                                           logger = NULL) {

  km <- nrow(Gamma_coeffs_matrix)
  V <- ncol(Gamma_coeffs_matrix)

  if (km != k_conditions * m_manifold_dim) {
    stop("Gamma_coeffs_matrix has incorrect number of rows")
  }

  Xi_raw <- matrix(0, m_manifold_dim, V)
  Beta_raw <- matrix(0, k_conditions, V)

  quality_metrics <- list(
    condition_numbers = numeric(V),
    svd_method = character(V),
    regularization_applied = logical(V),
    singular_value_gaps = numeric(V)
  )

  for (v in 1:V) {
    gamma_v <- Gamma_coeffs_matrix[, v]
    Gamma_mat <- matrix(gamma_v, nrow = k_conditions, ncol = m_manifold_dim, byrow = TRUE)

    if (all(abs(gamma_v) < .Machine$double.eps)) {
      Xi_raw[, v] <- 0
      Beta_raw[, v] <- 0
      quality_metrics$svd_method[v] <- "zero"
      next
    }

    gamma_scale <- max(abs(Gamma_mat))
    if (gamma_scale > 0) {
      Gamma_scaled <- Gamma_mat / gamma_scale
      cn <- kappa(Gamma_scaled)
      quality_metrics$condition_numbers[v] <- cn
      if (cn > max_condition_number) {
        reg_amount <- (cn / max_condition_number) * regularization_factor * .Machine$double.eps
        diag(Gamma_mat) <- diag(Gamma_mat) + reg_amount
        quality_metrics$regularization_applied[v] <- TRUE
        warning(sprintf("Voxel %d: Applied regularization due to condition number %.2e", v, cn))
      }
    }

    svd_result <- tryCatch({
      if (use_randomized_svd && k_conditions > 10 && m_manifold_dim > 10) {
        quality_metrics$svd_method[v] <- "randomized"
        if (requireNamespace("rsvd", quietly = TRUE)) {
          rsvd::rsvd(Gamma_mat, k = min(k_conditions, m_manifold_dim))
        } else {
          svd(Gamma_mat)
        }
      } else {
        quality_metrics$svd_method[v] <- "standard"
        svd(Gamma_mat)
      }
    }, error = function(e) {
      warning(sprintf("SVD failed for voxel %d: %s. Using fallback.", v, e$message))
      quality_metrics$svd_method[v] <- "fallback"
      if (all(is.finite(Gamma_mat)) && sum(Gamma_mat^2) > 0) {
        list(
          u = matrix(1/sqrt(k_conditions), k_conditions, 1),
          v = matrix(1/sqrt(m_manifold_dim), m_manifold_dim, 1),
          d = sqrt(sum(Gamma_mat^2))
        )
      } else {
        list(
          u = matrix(0, k_conditions, 1),
          v = matrix(0, m_manifold_dim, 1),
          d = 0
        )
      }
    })

    d <- svd_result$d

    if (length(d) > 1) {
      gaps <- diff(d) / d[-length(d)]
      quality_metrics$singular_value_gaps[v] <- max(abs(gaps))
      weights <- d / (d[1] + .Machine$double.eps)
      weights[weights < 0.01] <- 0
    } else {
      weights <- 1
    }

    if (length(d) > 0 && d[1] > .Machine$double.eps) {
      Xi_raw[, v] <- svd_result$v[, 1] * sqrt(d[1]) * weights[1]
      Beta_raw[, v] <- svd_result$u[, 1] * sqrt(d[1]) * weights[1]
    } else {
      Xi_raw[, v] <- 0
      Beta_raw[, v] <- 0
      quality_metrics$svd_method[v] <- "degenerate"
    }
  }

  n_regularized <- sum(quality_metrics$regularization_applied)
  if (n_regularized > 0) {
    msg <- sprintf("Applied regularization to %d/%d voxels", n_regularized, V)
    if (!is.null(logger)) logger$add(msg) else message(msg)
  }

  n_fallback <- sum(quality_metrics$svd_method == "fallback")
  if (n_fallback > 0) {
    msg <- sprintf("Used fallback SVD for %d/%d voxels", n_fallback, V)
    if (!is.null(logger)) logger$add(msg) else message(msg)
  }

  list(
    Xi_raw_matrix = Xi_raw,
    Beta_raw_matrix = Beta_raw,
    quality_metrics = quality_metrics
  )
}

#' Apply intrinsic identifiability constraints
#'
#' @param Xi_raw_matrix m x V matrix of raw manifold coordinates
#' @param Beta_raw_matrix k x V matrix of raw condition amplitudes
#' @param B_reconstructor_matrix p x m manifold reconstructor
#' @param h_ref_shape_vector p-length canonical HRF shape
#' @param ident_scale_method one of "l2_norm", "max_abs_val", "none"
#' @param zero_tol numeric tolerance for treating a reconstructed HRF as zero.
#'   Voxels with L2 norm or maximum absolute value below this threshold are
#'   zeroed in both \code{Xi_ident_matrix} and \code{Beta_ident_matrix}.
#' @param ident_sign_method Sign alignment method. Only
#'   "canonical_correlation" is currently supported.
#' @param consistency_check Logical; if TRUE, reprojection is used to verify
#'   alignment with the canonical HRF and voxels are flipped if necessary.
#' @return list with Xi_ident_matrix and Beta_ident_matrix
#' @export
apply_intrinsic_identifiability_core <- function(Xi_raw_matrix,
                                                 Beta_raw_matrix,
                                                 B_reconstructor_matrix,
                                                 h_ref_shape_vector,
                                                 ident_scale_method = c("l2_norm", "max_abs_val", "none"),
                                                 ident_sign_method = c("first_component", "canonical_correlation", "data_fit_correlation"),
                                                 zero_tol = 1e-8,
                                                 Y_proj_matrix = NULL,
                                                 X_condition_list_proj_matrices = NULL,
                                                 consistency_check = FALSE) {

  ident_scale_method <- match.arg(ident_scale_method)
  ident_sign_method <- match.arg(ident_sign_method)

  message(sprintf("Using '%s' for sign alignment", ident_sign_method))

  if (ident_sign_method == "data_fit_correlation") {
    if (is.null(Y_proj_matrix) || is.null(X_condition_list_proj_matrices)) {
      stop("data_fit_correlation method requires Y_proj_matrix and X_condition_list_proj_matrices")
    }
  }

  xi_ref_coord <- MASS::ginv(B_reconstructor_matrix) %*% h_ref_shape_vector

  V <- ncol(Xi_raw_matrix)
  Xi_ident <- matrix(0, nrow(Xi_raw_matrix), V)
  Beta_ident <- matrix(0, nrow(Beta_raw_matrix), V)

  for (v in seq_len(V)) {
    xi_v <- Xi_raw_matrix[, v]
    beta_v <- Beta_raw_matrix[, v]

    if (all(xi_v == 0)) {
      Xi_ident[, v] <- 0
      Beta_ident[, v] <- 0
      next
    }

    if (ident_sign_method == "canonical_correlation" || ident_sign_method == "first_component") {
      h_tmp <- B_reconstructor_matrix %*% xi_v
      
      # Verify dimensions match
      if (length(h_tmp) != length(h_ref_shape_vector)) {
        stop(sprintf("Dimension mismatch at voxel %d: reconstructed HRF has length %d but reference has length %d", 
                     v, length(h_tmp), length(h_ref_shape_vector)))
      }
      
      corr_ref <- tryCatch(cor(as.vector(h_tmp), as.vector(h_ref_shape_vector)), 
                          warning = function(w) NA, 
                          error = function(e) NA)
      if (is.na(corr_ref)) corr_ref <- 0
      if (abs(corr_ref) < 1e-3) {
        warning(sprintf("Voxel %d: canonical correlation near zero (%.3f)", v, corr_ref))
      }
      sgn <- sign(corr_ref)
      if (sgn == 0) sgn <- 1
      if (abs(corr_ref) < 1e-3 && !is.null(Y_proj_matrix) &&
          !is.null(X_condition_list_proj_matrices)) {
        best_r2 <- -Inf
        best_sgn <- 1
        for (sg in c(1, -1)) {
          xi_tmp <- xi_v * sg
          beta_tmp <- beta_v * sg
          h_tmp2 <- B_reconstructor_matrix %*% xi_tmp
          k2 <- length(X_condition_list_proj_matrices)
          X_design <- matrix(0, nrow(Y_proj_matrix), k2)
          for (c in 1:k2) {
            X_design[, c] <- X_condition_list_proj_matrices[[c]] %*% h_tmp2
          }
          y_pred <- X_design %*% beta_tmp
          y_true <- Y_proj_matrix[, v]
          r2 <- tryCatch(cor(as.vector(y_pred), as.vector(y_true))^2, 
                        warning = function(w) NA, 
                        error = function(e) NA)
          if (!is.na(r2) && r2 > best_r2) {
            best_r2 <- r2
            best_sgn <- sg
          }
        }
        if (best_r2 > 0) {
          sgn <- best_sgn
        } else {
          sgn <- sign(sum(h_tmp))
          if (sgn == 0) sgn <- 1
        }
      }
    } else if (ident_sign_method == "data_fit_correlation") {
      best_r2 <- -Inf
      best_sgn <- 1
      for (sg in c(1, -1)) {
        xi_tmp <- xi_v * sg
        beta_tmp <- beta_v * sg
        h_tmp <- B_reconstructor_matrix %*% xi_tmp
        k2 <- length(X_condition_list_proj_matrices)
        X_design <- matrix(0, nrow(Y_proj_matrix), k2)
        for (c in 1:k2) {
          X_design[, c] <- X_condition_list_proj_matrices[[c]] %*% h_tmp
        }
        y_pred <- X_design %*% beta_tmp
        y_true <- Y_proj_matrix[, v]
        r2 <- tryCatch(cor(as.vector(y_pred), as.vector(y_true))^2, 
                      warning = function(w) NA, 
                      error = function(e) NA)
        if (!is.na(r2) && r2 > best_r2) {
          best_r2 <- r2
          best_sgn <- sg
        }
      }
      sgn <- best_sgn
      if (!(best_r2 > 0)) {
        h_tmp <- B_reconstructor_matrix %*% xi_v
        corr_ref <- tryCatch(cor(as.vector(h_tmp), as.vector(h_ref_shape_vector)), 
                            warning = function(w) NA, 
                            error = function(e) NA)
        if (!is.na(corr_ref) && abs(corr_ref) >= 1e-3) {
          sgn <- sign(corr_ref)
        } else {
          sgn <- sign(sum(h_tmp))
          if (sgn == 0) sgn <- 1
        }
      }
    }

    xi_v <- xi_v * sgn
    beta_v <- beta_v * sgn

    hrf_v <- B_reconstructor_matrix %*% xi_v
    scale_val <- 1
    if (ident_scale_method == "l2_norm") {
      l2_norm <- sqrt(sum(hrf_v^2))
      if (l2_norm < zero_tol) {
        Xi_ident[, v] <- 0
        Beta_ident[, v] <- 0
        next
      }
      scale_val <- 1 / pmax(l2_norm, .Machine$double.eps)
    } else if (ident_scale_method == "max_abs_val") {
      max_abs <- max(abs(hrf_v))
      if (max_abs < zero_tol) {
        Xi_ident[, v] <- 0
        Beta_ident[, v] <- 0
        next
      }
      scale_val <- 1 / pmax(max_abs, .Machine$double.eps)
    }
    xi_out <- xi_v * scale_val
    beta_out <- beta_v / scale_val

    if (consistency_check) {
      hr_check <- B_reconstructor_matrix %*% xi_out
      corr_check <- tryCatch(cor(as.vector(hr_check), as.vector(h_ref_shape_vector)), 
                            warning = function(w) NA, 
                            error = function(e) NA)
      if (!is.na(corr_check) && corr_check < 0) {
        xi_out <- -xi_out
        beta_out <- -beta_out
      }
    }

    Xi_ident[, v] <- xi_out
    Beta_ident[, v] <- beta_out
  }

  list(Xi_ident_matrix = Xi_ident, Beta_ident_matrix = Beta_ident)
}


#' Construct voxel graph Laplacian
#'
#' @param voxel_coords_matrix V x 3 matrix of voxel coordinates
#' @param num_neighbors_Lsp number of nearest neighbours
#' @return sparse V x V graph Laplacian matrix
#' @export
make_voxel_graph_laplacian_core <- function(voxel_coords_matrix, num_neighbors_Lsp = 6,
                                            distance_engine = c("euclidean", "ann_euclidean"),
                                            ann_threshold = 10000) {
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
    distance_engine <- match.arg(distance_engine)
    if (distance_engine == "ann_euclidean" ||
        (distance_engine == "euclidean" && n_voxels > ann_threshold &&
         requireNamespace("RcppHNSW", quietly = TRUE))) {
      ann_res <- RcppHNSW::hnsw_knn(voxel_coords_matrix, k = num_neighbors_Lsp + 1)
      idx_mat <- ann_res$idx[, -1, drop = FALSE]
    } else {
      res <- knn_search_cpp(t(voxel_coords_matrix), t(voxel_coords_matrix), num_neighbors_Lsp + 1)
      idx_mat <- t(res$idx)[, -1, drop = FALSE]
    }
    i_idx <- rep(seq_len(n_voxels), each = ncol(idx_mat))
    j_idx <- as.vector(idx_mat)
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

  # Solve for all dimensions at once
  Xi_s_t <- Matrix::solve(A, t(Xi_ident_matrix))

  # Return regular matrix in m x V order
  t(as.matrix(Xi_s_t))
}



