#' Fast Least Squares Separate implementation
#'
#' Provides a straightforward LSS estimator that projects out fixed regressors
#' and then computes trial-wise betas using a numerically stable formula.
#'
#' @param Y n x V numeric matrix of BOLD data
#' @param dmat_base n x q matrix of baseline regressors (e.g., intercept, drift)
#' @param dmat_ran n x T matrix of trial regressors
#' @param dmat_fixed optional n x q2 matrix of additional fixed regressors
#' @return T x V matrix of trial-wise beta estimates
#' @keywords internal
lss_fast <- function(Y, dmat_base, dmat_ran, dmat_fixed = NULL) {
  if (!is.matrix(Y)) stop("Y must be a matrix")
  n_timepoints <- nrow(Y)
  if (nrow(dmat_base) != n_timepoints)
    stop("dmat_base dimension mismatch")
  if (nrow(dmat_ran) != n_timepoints)
    stop("dmat_ran dimension mismatch")
  if (!is.null(dmat_fixed) && nrow(dmat_fixed) != n_timepoints)
    stop("dmat_fixed dimension mismatch")

  X_base_fixed <- if (is.null(dmat_fixed)) dmat_base else cbind(dmat_base, dmat_fixed)
  P_base_fixed <- MASS::ginv(X_base_fixed)
  Q_base_fixed <- diag(n_timepoints) - X_base_fixed %*% P_base_fixed
  residual_data <- Q_base_fixed %*% Y
  Q_dmat_ran <- Q_base_fixed %*% dmat_ran
  lss_compute_r(Q_dmat_ran, residual_data)
}

#' Internal LSS computation used by \code{lss_fast}
#'
#' @param Q_dmat_ran n x T matrix of projected trial regressors
#' @param residual_data n x V matrix of data after projection
#' @return T x V matrix of beta estimates
#' @keywords internal
lss_compute_r <- function(Q_dmat_ran, residual_data) {
  n_events <- ncol(Q_dmat_ran)
  n_voxels <- ncol(residual_data)

  beta_matrix <- matrix(0, n_events, n_voxels)

  for (i in seq_len(n_events)) {
    c <- Q_dmat_ran[, i, drop = FALSE]
    c_norm <- drop(crossprod(c))
    if (c_norm < 1e-10) {
      beta_matrix[i, ] <- 0
      next
    }

    if (n_events == 1) {
      s <- c / c_norm
      beta_matrix[i, ] <- drop(crossprod(s, residual_data))
      next
    }

    b <- Q_dmat_ran[, -i, drop = FALSE]
    if (ncol(b) == 0 || all(colSums(b^2) < 1e-10)) {
      s <- c / c_norm
      beta_matrix[i, ] <- drop(crossprod(s, residual_data))
      next
    }

    b_sum <- rowSums(b)
    b <- matrix(b_sum, ncol = 1)
    b_norm <- drop(crossprod(b))
    if (b_norm < 1e-10) {
      s <- c / c_norm
      beta_matrix[i, ] <- drop(crossprod(s, residual_data))
      next
    }

    bc <- drop(crossprod(b, c))
    v <- c - (bc / b_norm) * b
    cvdot <- drop(crossprod(c, v))
    v_norm <- drop(crossprod(v))
    if (abs(cvdot) < 1e-8 * sqrt(c_norm * v_norm)) {
      s <- c / c_norm
    } else {
      s <- v / cvdot
    }
    beta_matrix[i, ] <- drop(crossprod(s, residual_data))
  }
  beta_matrix
}
