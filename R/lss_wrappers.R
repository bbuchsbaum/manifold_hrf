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
  hrfals::lss_mode_a(A = design_matrix, Y = data_matrix, C = trial_matrix, ...)
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
  hrfals::lss_mode_b(design_matrix, data_matrix, ...)
}
