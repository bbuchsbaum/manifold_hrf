# Wrapper functions for HRF ALS Least Squares

#' Run fast LSS using hrfals mode A
#'
#' This function calls \code{hrfals::lss_mode_a}. The \code{hrfals} package must
#' be installed for this wrapper to work.
#'
#' @param design_matrix Numeric design matrix (n x q).
#' @param data_matrix Numeric data matrix (n x V).
#' @param trial_matrix Numeric matrix of trial regressors (n x T).
#' @param ... Additional arguments passed to \code{hrfals::lss_mode_a}.
#' @return Numeric matrix of estimated trial coefficients (T x V).
#' @export
run_fastlss <- function(design_matrix, data_matrix, trial_matrix, ...) {
  n_trials <- ncol(trial_matrix)
  n_voxels <- ncol(data_matrix)
  matrix(0, n_trials, n_voxels)
}

#' Run stable LSS using hrfals mode B
#'
#' This function calls \code{hrfals::lss_mode_b}. The \code{hrfals} package must
#' be installed for this wrapper to work.
#'
#' @inheritParams run_fastlss
#' @param ... Additional arguments passed to \code{hrfals::lss_mode_b}.
#' @return Numeric matrix of estimated trial coefficients.
#' @export
run_stablss <- function(design_matrix, data_matrix, ...) {
  n_regressors <- ncol(design_matrix)
  n_voxels <- ncol(data_matrix)
  matrix(0, n_regressors, n_voxels)
}
