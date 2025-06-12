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

#' Create FIR (Finite Impulse Response) Basis
#'
#' Creates a finite impulse response basis set for HRF estimation.
#' This is a more standard alternative to HRF_RAW_EVENT_BASIS that
#' doesn't depend on fmrireg internals.
#'
#' @param n_basis Number of basis functions (time points)
#' @param TR Repetition time in seconds
#' @param name Optional name for the basis
#' @return An HRF object compatible with fmrireg
#' @export
create_fir_basis <- function(n_basis, TR, name = NULL) {
  # Input validation
  n_basis <- as.integer(n_basis)
  if (n_basis < 1) {
    stop("n_basis must be at least 1")
  }
  
  if (!is.numeric(TR) || TR <= 0) {
    stop("TR must be a positive number")
  }
  
  # Total duration of the HRF
  span <- n_basis * TR
  
  # Create the basis evaluation function
  basis_fun <- function(t) {
    if (length(t) == 0) {
      return(matrix(0, nrow = 0, ncol = n_basis))
    }
    
    # Initialize output matrix
    out <- matrix(0, nrow = length(t), ncol = n_basis)
    
    # For each time point, determine which basis function is active
    # FIR basis uses indicator functions for each time bin
    for (i in 1:n_basis) {
      # Time window for this basis function
      t_start <- (i - 1) * TR
      t_end <- i * TR
      
      # Set to 1 where t falls in this window
      in_window <- t >= t_start & t < t_end
      out[in_window, i] <- 1
    }
    
    # Handle the edge case for the last basis function
    # Include t == span in the last bin
    if (n_basis > 0) {
      out[t == span, n_basis] <- 1
    }
    
    return(out)
  }
  
  # Set name if not provided
  if (is.null(name)) {
    name <- sprintf("FIR_%d_TR%.2f", n_basis, TR)
  }
  
  # Create HRF object using fmrireg
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    stop("Package 'fmrireg' is required. Install it with: remotes::install_github('bbuchsbaum/fmrireg')")
  }
  
  hrf_obj <- fmrireg::as_hrf(
    basis_fun, 
    name = name, 
    nbasis = n_basis, 
    span = span,
    params = list(n_basis = n_basis, TR = TR)
  )
  
  return(hrf_obj)
}
