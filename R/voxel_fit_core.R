# Unified voxel fitting core API
# Combines standard and robust implementations

#' Project out confounds from data and design matrices
#'
#' @param Y_data_matrix n x V matrix of fMRI time series data
#' @param X_list_of_matrices List of n x p design matrices for each condition
#' @param Z_confounds_matrix Optional n x q matrix of confound regressors
#' @return List with Y_proj_matrix and X_list_proj_matrices
#' @export
project_out_confounds_core <- function(Y_data_matrix,
                                      X_list_of_matrices,
                                      Z_confounds_matrix = NULL) {
  if (!is.matrix(Y_data_matrix)) {
    stop("Y_data_matrix must be a matrix")
  }
  if (anyNA(Y_data_matrix)) {
    stop("Y_data_matrix must not contain NA values")
  }
  n <- nrow(Y_data_matrix)
  validate_design_matrix_list(X_list_of_matrices, n_timepoints = n)
  if (is.null(Z_confounds_matrix)) {
    return(list(
      Y_proj_matrix = Y_data_matrix,
      X_list_proj_matrices = X_list_of_matrices
    ))
  }
  validate_confounds_matrix(Z_confounds_matrix, n_timepoints = n)
  if (ncol(Z_confounds_matrix) >= n) {
    stop("Z_confounds_matrix has too many columns (must be less than number of timepoints)")
  }
  svd_Z <- svd(Z_confounds_matrix)
  tol_svd <- max(dim(Z_confounds_matrix)) * max(svd_Z$d) * .Machine$double.eps
  rank_Z <- sum(svd_Z$d > tol_svd)
  if (rank_Z < ncol(Z_confounds_matrix)) {
    warning("Z_confounds_matrix is rank deficient; using independent columns only")
  }
  Q_Z <- svd_Z$u[, seq_len(rank_Z), drop = FALSE]
  Y_proj_matrix <- Y_data_matrix - Q_Z %*% crossprod(Q_Z, Y_data_matrix)
  X_list_proj_matrices <- lapply(X_list_of_matrices, function(X) {
    X - Q_Z %*% crossprod(Q_Z, X)
  })
  list(Y_proj_matrix = Y_proj_matrix, X_list_proj_matrices = X_list_proj_matrices)
}

#' Transform designs to manifold basis
#'
#' @param X_condition_list_proj_matrices List of n x p projected design matrices
#' @param B_reconstructor_matrix p x m manifold basis matrix
#' @return List of n x m transformed design matrices
#' @export
transform_designs_to_manifold_basis_core <- function(X_condition_list_proj_matrices,
                                                    B_reconstructor_matrix) {
  if (!is.list(X_condition_list_proj_matrices)) {
    stop("X_condition_list_proj_matrices must be a list")
  }
  if (!is.matrix(B_reconstructor_matrix)) {
    stop("B_reconstructor_matrix must be a matrix")
  }
  if (length(X_condition_list_proj_matrices) == 0) {
    stop("X_condition_list_proj_matrices cannot be empty")
  }
  p <- nrow(B_reconstructor_matrix)
  m <- ncol(B_reconstructor_matrix)
  for (i in seq_along(X_condition_list_proj_matrices)) {
    if (!is.matrix(X_condition_list_proj_matrices[[i]])) {
      stop(sprintf("X_condition_list_proj_matrices[[%d]] must be a matrix", i))
    }
    if (ncol(X_condition_list_proj_matrices[[i]]) != p) {
      stop(sprintf(
        "X_condition_list_proj_matrices[[%d]] has %d columns but B_reconstructor_matrix has %d rows",
        i, ncol(X_condition_list_proj_matrices[[i]]), p
      ))
    }
  }
  lapply(X_condition_list_proj_matrices, function(X_i) X_i %*% B_reconstructor_matrix)
}

#' Solve GLM for gamma coefficients
#'
#' @param Z_list_of_matrices List of n x m transformed design matrices
#' @param Y_proj_matrix n x V projected data matrix
#' @param lambda_gamma Ridge regularization parameter
#' @param orthogonal_approx_flag Whether to use orthogonal approximation
#' @return (k*m) x V matrix of gamma coefficients
#' @export
solve_glm_for_gamma_core <- function(Z_list_of_matrices,
                                    Y_proj_matrix,
                                    lambda_gamma,
                                    orthogonal_approx_flag = FALSE) {
  if (!is.list(Z_list_of_matrices)) {
    stop("Z_list_of_matrices must be a list")
  }
  if (length(Z_list_of_matrices) == 0) {
    stop("Z_list_of_matrices cannot be empty")
  }
  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }
  lambda_gamma <- .validate_and_standardize_lambda(lambda_gamma, "lambda_gamma")
  if (!is.logical(orthogonal_approx_flag) || length(orthogonal_approx_flag) != 1) {
    stop("orthogonal_approx_flag must be a single logical value")
  }
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  k <- length(Z_list_of_matrices)
  if (!is.matrix(Z_list_of_matrices[[1]])) {
    stop("Z_list_of_matrices[[1]] must be a matrix")
  }
  m <- ncol(Z_list_of_matrices[[1]])
  for (i in seq_along(Z_list_of_matrices)) {
    if (!is.matrix(Z_list_of_matrices[[i]])) {
      stop(sprintf("Z_list_of_matrices[[%d]] must be a matrix", i))
    }
    if (nrow(Z_list_of_matrices[[i]]) != n) {
      stop(sprintf("Z_list_of_matrices[[%d]] must have %d rows to match Y_proj_matrix", i, n))
    }
    if (ncol(Z_list_of_matrices[[i]]) != m) {
      stop(sprintf("All Z matrices must have the same number of columns. Z[[%d]] has %d columns but Z[[1]] has %d",
                   i, ncol(Z_list_of_matrices[[i]]), m))
    }
  }
  X_tilde <- do.call(cbind, Z_list_of_matrices)
  XtX <- crossprod(X_tilde)
  if (orthogonal_approx_flag) {
    XtX_approx <- matrix(0, nrow = k * m, ncol = k * m)
    for (i in 1:k) {
      idx <- ((i - 1) * m + 1):(i * m)
      XtX_approx[idx, idx] <- XtX[idx, idx]
    }
    XtX <- XtX_approx
  }
  XtX_reg <- XtX + lambda_gamma * diag(k * m)
  
  # Check for perfect collinearity before attempting solve
  if (lambda_gamma == 0) {
    # Check rank of XtX for perfect collinearity
    svd_XtX <- svd(XtX, nu = 0, nv = 0)
    tol <- max(dim(XtX)) * max(svd_XtX$d) * .Machine$double.eps * 100
    rank_XtX <- sum(svd_XtX$d > tol)
    if (rank_XtX < k * m) {
      stop("Design matrix is singular (perfectly collinear columns)")
    }
  }
  
  cond_num <- kappa(XtX_reg, exact = FALSE)
  if (cond_num > 1e10) {
    warning(sprintf("XtX_reg has high condition number (%.2e). Consider increasing lambda_gamma.", cond_num))
  }
  XtY <- crossprod(X_tilde, Y_proj_matrix)
  
  # Try to solve, handling near-singular cases
  tryCatch({
    solve(XtX_reg, XtY, tol = .Machine$double.eps)
  }, error = function(e) {
    if (grepl("singular|computationally singular", e$message, ignore.case = TRUE) && lambda_gamma == 0) {
      # For near-singular cases with no regularization, 
      # use SVD-based pseudoinverse to get potentially large coefficients
      warning("Near-singular matrix detected. Using pseudoinverse.")
      svd_XtX <- svd(XtX_reg)
      # Keep only non-zero singular values
      tol <- max(dim(XtX_reg)) * max(svd_XtX$d) * .Machine$double.eps
      pos <- svd_XtX$d > tol
      if (sum(pos) == 0) {
        # Completely singular - return zeros
        return(matrix(0, nrow = k * m, ncol = V))
      }
      # Compute pseudoinverse
      d_inv <- numeric(length(svd_XtX$d))
      d_inv[pos] <- 1 / svd_XtX$d[pos]
      XtX_pinv <- svd_XtX$v %*% (d_inv * t(svd_XtX$u))
      XtX_pinv %*% XtY
    } else {
      stop(e)
    }
  })
}

# internal standard SVD implementation
.extract_xi_beta_raw_svd_standard <- function(Gamma_coeffs_matrix,
                                             m_manifold_dim,
                                             k_conditions) {
  V <- ncol(Gamma_coeffs_matrix)
  Xi_raw <- matrix(0, m_manifold_dim, V)
  Beta_raw <- matrix(0, k_conditions, V)
  for (v in seq_len(V)) {
    gamma_v <- Gamma_coeffs_matrix[, v]
    
    # Check for NA/NaN values
    if (any(!is.finite(gamma_v))) {
      warning(sprintf("Voxel %d: Contains NA/NaN/Inf values, skipping", v))
      next
    }
    
    Gv <- matrix(gamma_v, nrow = k_conditions,
                 ncol = m_manifold_dim, byrow = TRUE)
    
    # Check for zero matrix
    if (all(abs(Gv) < .Machine$double.eps)) {
      next
    }
    
    # Perform SVD with error handling
    sv <- tryCatch({
      svd(Gv)
    }, error = function(e) {
      warning(sprintf("Voxel %d: SVD failed - %s", v, e$message))
      return(NULL)
    })
    
    if (is.null(sv) || length(sv$d) == 0 || sv$d[1] < .Machine$double.eps) next
    
    Beta_raw[, v] <- sv$u[, 1] * sqrt(sv$d[1])
    Xi_raw[, v] <- sv$v[, 1] * sqrt(sv$d[1])
  }
  list(Xi_raw_matrix = Xi_raw, Beta_raw_matrix = Beta_raw)
}

# robust SVD implementation
.extract_xi_beta_raw_svd_robust_impl <- function(Gamma_coeffs_matrix,
                                                m_manifold_dim,
                                                k_conditions,
                                                regularization_factor = 10,
                                                max_condition_number = 1e8,
                                                use_randomized_svd = FALSE,
                                                verbose_warnings = FALSE,
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
    
    # Check for NA/NaN/Inf values
    if (any(!is.finite(gamma_v))) {
      Xi_raw[, v] <- 0
      Beta_raw[, v] <- 0
      quality_metrics$svd_method[v] <- "non-finite"
      if (verbose_warnings) {
        warning(sprintf("Voxel %d: Contains non-finite values, skipping", v))
      }
      next
    }
    
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
      svd_check <- svd(Gamma_mat, nu = 0, nv = 0)
      n_sig <- sum(svd_check$d > max(svd_check$d) * .Machine$double.eps * 100)
      if (cn > max_condition_number && n_sig > 1) {
        reg_amount <- (cn / max_condition_number) * regularization_factor * .Machine$double.eps
        diag(Gamma_mat) <- diag(Gamma_mat) + reg_amount
        quality_metrics$regularization_applied[v] <- TRUE
        if (verbose_warnings) {
          warning(sprintf("Voxel %d: Applied regularization due to condition number %.2e", v, cn))
        }
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
        list(u = matrix(1/sqrt(k_conditions), k_conditions, 1),
             v = matrix(1/sqrt(m_manifold_dim), m_manifold_dim, 1),
             d = sqrt(sum(Gamma_mat^2)))
      } else {
        list(u = matrix(0, k_conditions, 1),
             v = matrix(0, m_manifold_dim, 1),
             d = 0)
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
  list(Xi_raw_matrix = Xi_raw, Beta_raw_matrix = Beta_raw, quality_metrics = quality_metrics)
}

#' Extract Xi and Beta via SVD
#'
#' @param Gamma_coeffs_matrix (km) x V matrix of coefficients
#' @param m_manifold_dim Number of manifold dimensions
#' @param k_conditions Number of conditions
#' @param method Use 'standard' or 'robust' SVD
#' @param ... Additional arguments passed to the robust implementation
#' @export
extract_xi_beta_raw_svd_core <- function(Gamma_coeffs_matrix,
                                        m_manifold_dim,
                                        k_conditions,
                                        method = c("standard", "robust"),
                                        ...) {
  method <- match.arg(method)
  if (!is.matrix(Gamma_coeffs_matrix)) {
    stop("Gamma_coeffs_matrix must be a matrix")
  }
  
  # Validate m and k
  if (!is.numeric(m_manifold_dim) || length(m_manifold_dim) != 1 || 
      m_manifold_dim <= 0 || m_manifold_dim != round(m_manifold_dim)) {
    stop("m_manifold_dim must be a positive integer")
  }
  if (!is.numeric(k_conditions) || length(k_conditions) != 1 || 
      k_conditions <= 0 || k_conditions != round(k_conditions)) {
    stop("k_conditions must be a positive integer")
  }
  
  expected_rows <- k_conditions * m_manifold_dim
  if (nrow(Gamma_coeffs_matrix) != expected_rows) {
    stop(sprintf("Gamma_coeffs_matrix has %d rows but expected %d (k * m)",
                 nrow(Gamma_coeffs_matrix), expected_rows))
  }
  
  # Check for NA/NaN/Inf values in the matrix
  if (any(!is.finite(Gamma_coeffs_matrix))) {
    n_bad <- sum(!is.finite(Gamma_coeffs_matrix))
    n_total <- length(Gamma_coeffs_matrix)
    warning(sprintf("Gamma_coeffs_matrix contains %d non-finite values (%.1f%%). These will be handled per voxel.",
                    n_bad, 100 * n_bad / n_total))
  }
  
  if (method == "robust") {
    .extract_xi_beta_raw_svd_robust_impl(Gamma_coeffs_matrix, m_manifold_dim, k_conditions, ...)
  } else {
    .extract_xi_beta_raw_svd_standard(Gamma_coeffs_matrix, m_manifold_dim, k_conditions)
  }
}

#' Wrapper for robust SVD extraction
#' 
#' @param Gamma_coeffs_matrix (km) x V matrix of coefficients
#' @param m_manifold_dim Number of manifold dimensions
#' @param k_conditions Number of conditions
#' @param ... Additional arguments passed to the robust implementation
#' @return List with Xi_raw_matrix and Beta_raw_matrix
#' @export
extract_xi_beta_raw_svd_robust <- function(Gamma_coeffs_matrix,
                                          m_manifold_dim,
                                          k_conditions,
                                          ...) {
  extract_xi_beta_raw_svd_core(Gamma_coeffs_matrix, m_manifold_dim, k_conditions,
                               method = "robust", ...)
}

#' Apply intrinsic identifiability constraints
#'
#' @param Xi_raw_matrix m x V matrix of raw manifold coordinates
#' @param Beta_raw_matrix k x V matrix of raw condition amplitudes
#' @param B_reconstructor_matrix p x m manifold basis matrix
#' @param h_ref_shape_vector Reference HRF shape vector of length p
#' @param ident_scale_method Scaling method: "l2_norm", "max_abs_val", or "none"
#' @param ident_sign_method Sign alignment method: "canonical_correlation" or "data_fit_correlation"
#' @param zero_tol Tolerance for considering values as zero
#' @param correlation_threshold Threshold for correlation-based alignment
#' @param Y_proj_matrix Optional n x V projected data matrix (for data_fit_correlation)
#' @param X_condition_list_proj_matrices Optional list of projected design matrices
#' @param consistency_check Whether to perform consistency checks
#' @param n_jobs Number of parallel jobs
#' @param verbose Whether to print progress messages
#' @return List with Xi_ident_matrix and Beta_ident_matrix
#' @export
apply_intrinsic_identifiability_core <- function(Xi_raw_matrix,
                                                Beta_raw_matrix,
                                                B_reconstructor_matrix,
                                                h_ref_shape_vector,
                                                ident_scale_method = "l2_norm",
                                                ident_sign_method = "canonical_correlation",
                                                zero_tol = 1e-8,
                                                correlation_threshold = 1e-3,
                                                Y_proj_matrix = NULL,
                                                X_condition_list_proj_matrices = NULL,
                                                consistency_check = FALSE,
                                                n_jobs = 1,
                                                verbose = FALSE) {
  if (!is.matrix(Xi_raw_matrix)) stop("Xi_raw_matrix must be a matrix")
  if (!is.matrix(Beta_raw_matrix)) stop("Beta_raw_matrix must be a matrix")
  if (!is.matrix(B_reconstructor_matrix)) stop("B_reconstructor_matrix must be a matrix")
  if (!is.numeric(h_ref_shape_vector) || !is.vector(h_ref_shape_vector)) {
    stop("h_ref_shape_vector must be a numeric vector")
  }
  m <- nrow(Xi_raw_matrix)
  V <- ncol(Xi_raw_matrix)
  k <- nrow(Beta_raw_matrix)
  p <- nrow(B_reconstructor_matrix)
  if (ncol(Beta_raw_matrix) != V) {
    stop("Xi_raw_matrix and Beta_raw_matrix must have the same number of columns")
  }
  if (ncol(B_reconstructor_matrix) != m) {
    stop("B_reconstructor_matrix must have m columns to match Xi dimension")
  }
  if (length(h_ref_shape_vector) != p) {
    stop("h_ref_shape_vector must have length p to match B_reconstructor rows")
  }
  valid_scale_methods <- c("l2_norm", "max_abs_val", "none")
  if (!ident_scale_method %in% valid_scale_methods) {
    stop("ident_scale_method must be one of: ", paste(valid_scale_methods, collapse = ", "))
  }
  valid_sign_methods <- c("canonical_correlation", "data_fit_correlation")
  if (!ident_sign_method %in% valid_sign_methods) {
    stop("ident_sign_method must be one of: ", paste(valid_sign_methods, collapse = ", "))
  }
  message(sprintf("Using '%s' for sign alignment", ident_sign_method))
  if (ident_sign_method == "canonical_correlation") {
    result <- apply_identifiability_vectorized(
      Xi_raw = Xi_raw_matrix,
      Beta_raw = Beta_raw_matrix,
      B = B_reconstructor_matrix,
      h_ref = h_ref_shape_vector,
      h_mode = "max_correlation",
      scale_method = ident_scale_method
    )
    return(result)
  } else {
    if (is.null(Y_proj_matrix) || is.null(X_condition_list_proj_matrices)) {
      stop("data_fit_correlation method requires Y_proj_matrix and X_condition_list_proj_matrices")
    }
    config <- list(
      sign_method = ident_sign_method,
      scale_method = ident_scale_method,
      zero_tol = zero_tol,
      correlation_threshold = correlation_threshold,
      consistency_check = consistency_check
    )
    voxel_fun <- function(vx) {
      .process_voxel_identifiability(
        vx = vx,
        Xi_raw = Xi_raw_matrix,
        Beta_raw = Beta_raw_matrix,
        B_reconstructor = B_reconstructor_matrix,
        h_ref = h_ref_shape_vector,
        config = config,
        Y_proj = Y_proj_matrix,
        X_list = X_condition_list_proj_matrices,
        verbose = verbose
      )
    }
    res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
    Xi_ident_matrix <- do.call(cbind, lapply(res_list, "[[", "xi"))
    Beta_ident_matrix <- do.call(cbind, lapply(res_list, "[[", "beta"))
    list(Xi_ident_matrix = Xi_ident_matrix, Beta_ident_matrix = Beta_ident_matrix)
  }
}

#' Smart initialization for voxelwise fit
#'
#' @param Y_data n x V matrix of fMRI data
#' @param X_condition_list List of n x p design matrices for each condition
#' @param hrf_canonical Canonical HRF shape vector of length p
#' @param use_spatial_clusters Whether to use spatial clustering for initialization
#' @param voxel_coords V x 3 matrix of voxel coordinates (if use_spatial_clusters = TRUE)
#' @param m_manifold_dim Target manifold dimension
#' @return List with initialization values
#' @export
smart_initialize <- function(Y_data, X_condition_list, hrf_canonical,
                           use_spatial_clusters = TRUE,
                           voxel_coords = NULL,
                           m_manifold_dim = 5) {
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  k <- length(X_condition_list)
  message("Computing smart initialization...")
  
  # More memory-efficient approach: compute X_canonical one column at a time
  # and immediately compute XtX without storing full X_canonical
  XtX <- matrix(0, k, k)
  XtY <- matrix(0, k, V)
  
  # Process conditions one at a time to reduce memory usage
  for (c in 1:k) {
    X_c <- as.vector(X_condition_list[[c]] %*% hrf_canonical)
    
    # Update XtX
    for (c2 in c:k) {
      X_c2 <- if (c2 == c) X_c else as.vector(X_condition_list[[c2]] %*% hrf_canonical)
      XtX[c, c2] <- XtX[c2, c] <- crossprod(X_c, X_c2)
    }
    
    # Update XtY
    XtY[c, ] <- crossprod(X_c, Y_data)
  }
  
  # Regularize and invert
  XtX_reg <- XtX + diag(0.01, k)
  XtX_inv <- solve(XtX_reg)
  
  # Compute Beta_init directly from XtY
  Beta_init <- XtX_inv %*% XtY
  
  # Compute R2 in chunks to reduce memory usage
  R2_init <- numeric(V)
  chunk_size <- min(1000, V)
  n_chunks <- ceiling(V / chunk_size)
  
  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, V)
    chunk_voxels <- start_idx:end_idx
    
    # Reconstruct X_canonical for this chunk only
    X_canonical_chunk <- matrix(0, n, k)
    for (c in 1:k) {
      X_canonical_chunk[, c] <- X_condition_list[[c]] %*% hrf_canonical
    }
    
    for (v in chunk_voxels) {
      y_pred <- X_canonical_chunk %*% Beta_init[, v]
      y_v <- Y_data[, v]
      ss_tot <- sum((y_v - mean(y_v))^2)
      ss_res <- sum((y_v - y_pred)^2)
      
      if (ss_tot < .Machine$double.eps) {
        R2_init[v] <- 0
      } else {
        R2_init[v] <- 1 - ss_res / ss_tot
        R2_init[v] <- max(0, min(1, R2_init[v]))
      }
    }
  }
  good_voxels <- which(R2_init > quantile(R2_init, 0.75, na.rm = TRUE))
  if (use_spatial_clusters && !is.null(voxel_coords) && length(good_voxels) > 10) {
    n_clusters <- min(20, length(good_voxels) / 5)
    kmeans_result <- kmeans(voxel_coords[good_voxels, ], centers = n_clusters)
    cluster_centers <- kmeans_result$centers
    # More memory-efficient nearest cluster assignment
    nearest_cluster <- integer(V)
    chunk_size <- min(1000, V)
    n_chunks <- ceiling(V / chunk_size)
    
    for (chunk in 1:n_chunks) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, V)
      chunk_voxels <- start_idx:end_idx
      
      # Process chunk of voxels
      for (v in chunk_voxels) {
        distances <- rowSums((cluster_centers - matrix(voxel_coords[v, ], n_clusters, 3, byrow = TRUE))^2)
        nearest_cluster[v] <- which.min(distances)
      }
    }
    message(sprintf("Using %d spatial clusters for initialization", n_clusters))
  } else {
    nearest_cluster <- rep(1, V)
  }
  Xi_init <- matrix(rnorm(m_manifold_dim * V, sd = 0.1), m_manifold_dim, V)
  beta_scale <- apply(Beta_init, 2, function(x) sqrt(sum(x^2)))
  Xi_init <- Xi_init * rep(beta_scale, each = m_manifold_dim)
  list(
    Xi_init = Xi_init,
    Beta_init = Beta_init,
    R2_init = R2_init,
    good_voxels = good_voxels,
    nearest_cluster = nearest_cluster
  )
}

#' Detect and reinitialize stuck voxels
#'
#' @keywords internal
fix_stuck_voxels <- function(Xi_current, Xi_previous,
                           variance_threshold = 1e-6,
                           reinit_sd = 0.1) {
  m <- nrow(Xi_current)
  V <- ncol(Xi_current)
  if (!is.null(Xi_previous)) {
    Xi_change <- Xi_current - Xi_previous
    change_variance <- apply(Xi_change, 2, var)
    stuck <- which(change_variance < variance_threshold)
    if (length(stuck) > 0) {
      message(sprintf("Reinitializing %d stuck voxels", length(stuck)))
      Xi_current[, stuck] <- matrix(rnorm(m * length(stuck), sd = reinit_sd),
                                   m, length(stuck))
    }
  }
  Xi_current
}

#' Construct voxel graph Laplacian
#'
#' @param voxel_coords_matrix V x 3 matrix of voxel coordinates
#' @param num_neighbors_Lsp Number of nearest neighbors for graph construction
#' @param distance_engine Distance computation method: "euclidean" or "ann_euclidean"
#' @param ann_threshold Threshold for switching to approximate nearest neighbors
#' @return Sparse Laplacian matrix
#' @export
make_voxel_graph_laplacian_core <- function(voxel_coords_matrix, num_neighbors_Lsp = 6,
                                            distance_engine = c("euclidean", "ann_euclidean"),
                                            ann_threshold = 10000) {
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
  if (n_voxels <= num_neighbors_Lsp) {
    warning(sprintf("Requested %d neighbors but only %d other voxels available. Creating fully connected graph.",
                    num_neighbors_Lsp, n_voxels - 1))
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
