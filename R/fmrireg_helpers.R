# fmrireg helper functions for manifoldhrf integration

#' Create an HRF basis that encodes raw event time courses
#'
#' @description
#' This function produces an `fmrireg::HRF` object representing
#' a series of delta functions (or narrow boxcars) used to
#' extract raw event time courses of length `p_length`.
#' Each basis function corresponds to a single sample offset
#' from the event onset.
#'
#' @param p_length Integer number of samples in the raw HRF.
#' @param TR_sample Numeric sample spacing of the HRF in seconds.
#' @param name Optional name of the HRF object.
#'
#' @return An object of class `HRF` from \pkg{fmrireg}.
#' @keywords internal
HRF_RAW_EVENT_BASIS <- function(p_length, TR_sample, name = NULL) {
  nbasis <- as.integer(p_length)
  span <- nbasis * TR_sample

  basis_fun <- function(t) {
    if (length(t) == 0) {
      return(matrix(0, nrow = 0, ncol = nbasis))
    }
    out <- matrix(0, nrow = length(t), ncol = nbasis)
    idx <- floor(t / TR_sample) + 1
    valid <- which(t >= 0 & t < span & idx >= 1 & idx <= nbasis)
    if (length(valid) > 0) {
      out[cbind(valid, idx[valid])] <- 1
    }
    out
  }

  hrf_name <- name %||% sprintf("RawEventBasis_p%d_TR%.2f", nbasis, TR_sample)
  fmrireg::as_hrf(basis_fun, name = hrf_name, nbasis = nbasis, span = span,
                  params = list(p_length = p_length, TR_sample = TR_sample))
}

#' Create a token for a factor level following fmrireg conventions
#'
#' @param factor_name Name of the factor variable
#' @param level Level value
#'
#' @return Character string combining the factor name and level
#' @keywords internal
private_level_token <- function(factor_name, level) {
  paste0(make.names(factor_name), ".", make.names(level))
}
