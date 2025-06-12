# Manifold and HRF Utility Functions

# Soundness Improvements for M-HRF-LSS
# Implementation of SOUND-* improvements

#' Check HRF Library Quality
#'
#' Evaluates the quality and diversity of an HRF library matrix
#'
#' @param L_library_matrix p x N matrix of HRF shapes
#' @return List with quality metrics and flags
#' @keywords internal
check_hrf_library_quality <- function(L_library_matrix) {
  
  p <- nrow(L_library_matrix)
  N <- ncol(L_library_matrix)
  
  quality <- list()
  
  # Check for duplicate HRFs
  # Handle case where all HRFs are constant
  cor_result <- tryCatch({
    cor_matrix <- cor(L_library_matrix)
    diag(cor_matrix) <- 0
    list(
      max_cor = max(abs(cor_matrix)),
      n_near_duplicates = sum(abs(cor_matrix) > 0.99) / 2
    )
  }, error = function(e) {
    # If correlation fails (e.g., zero variance), treat as degenerate
    list(
      max_cor = 0,
      n_near_duplicates = 0
    )
  })
  
  max_cor <- cor_result$max_cor
  n_near_duplicates <- cor_result$n_near_duplicates
  
  quality$has_duplicates <- n_near_duplicates > 0
  quality$n_duplicates <- n_near_duplicates
  quality$max_correlation <- max_cor
  
  # Check condition number
  svd_L <- svd(L_library_matrix)
  condition_number <- svd_L$d[1] / svd_L$d[min(p, N)]
  quality$condition_number <- condition_number
  quality$is_ill_conditioned <- condition_number > 1e10
  
  # Check for all-zero or constant HRFs
  zero_hrfs <- apply(L_library_matrix, 2, function(x) all(x == 0))
  constant_hrfs <- apply(L_library_matrix, 2, function(x) sd(x) < .Machine$double.eps)
  
  quality$n_zero_hrfs <- sum(zero_hrfs)
  quality$n_constant_hrfs <- sum(constant_hrfs)
  quality$has_degenerate_hrfs <- any(zero_hrfs | constant_hrfs)
  
  # Check diversity (mean pairwise distance)
  if (N > 1) {
    distances <- dist(t(L_library_matrix))
    quality$mean_diversity <- mean(distances)
    quality$min_diversity <- min(distances)
    quality$is_low_diversity <- quality$min_diversity < 0.01
  } else {
    quality$mean_diversity <- NA
    quality$min_diversity <- NA
    quality$is_low_diversity <- TRUE
  }
  
  # Overall quality flag
  quality$is_good_quality <- !quality$has_duplicates && 
                            !quality$is_ill_conditioned && 
                            !quality$has_degenerate_hrfs && 
                            !quality$is_low_diversity
  
  return(quality)
}


#' Remove Duplicate HRFs from Library
#'
#' Removes near-duplicate HRFs based on correlation threshold
#'
#' @param L_library_matrix p x N matrix of HRF shapes
#' @param cor_threshold Correlation threshold for duplicates (default 0.99)
#' @return Cleaned matrix with duplicates removed
#' @keywords internal
remove_duplicate_hrfs <- function(L_library_matrix, cor_threshold = 0.99) {
  
  N <- ncol(L_library_matrix)
  if (N <= 1) return(L_library_matrix)
  
  # Compute correlations
  cor_matrix <- cor(L_library_matrix)
  
  # Find which HRFs to keep
  keep <- rep(TRUE, N)
  
  for (i in 1:(N-1)) {
    if (keep[i]) {
      # Mark duplicates of HRF i for removal
      duplicates <- which(cor_matrix[i, (i+1):N] > cor_threshold) + i
      if (length(duplicates) > 0) {
        keep[duplicates] <- FALSE
      }
    }
  }
  
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message(sprintf("Removed %d duplicate HRFs (correlation > %.2f)", 
                   n_removed, cor_threshold))
  }
  
  return(L_library_matrix[, keep, drop = FALSE])
}


#' Compute PCA Basis as Fallback
#'
#' Computes PCA-based reconstructor when manifold construction fails
#'
#' @param L_library_matrix p x N matrix of HRF shapes
#' @param m_target Target dimensionality
#' @param min_variance Minimum variance to retain
#' @return List with B_reconstructor_matrix and metadata
#' @keywords internal
compute_pca_fallback <- function(L_library_matrix, m_target, min_variance = 0.95) {
  
  message("Using PCA fallback for manifold construction")
  
  p <- nrow(L_library_matrix)
  N <- ncol(L_library_matrix)
  
  # Center the HRFs
  L_centered <- L_library_matrix - rowMeans(L_library_matrix)
  
  # Compute SVD
  svd_result <- svd(L_centered)
  
  # Determine dimensionality based on variance
  var_explained <- svd_result$d^2 / sum(svd_result$d^2)
  cum_var <- cumsum(var_explained)
  m_auto <- which(cum_var >= min_variance)[1]
  
  # Use minimum of target and auto-selected
  m_final <- min(m_target, m_auto, length(svd_result$d))
  
  # Ensure m_final is valid
  if (is.na(m_final) || m_final < 1) {
    warning("No valid components found in PCA fallback")
    m_final <- 1
  }
  
  # Ensure we explain at least min_variance
  if (m_final > 0 && m_final <= length(cum_var) && 
      !is.na(cum_var[m_final]) && cum_var[m_final] < min_variance && 
      m_final < length(svd_result$d)) {
    m_final <- m_auto
  }
  
  # Extract components
  if (m_final > 0 && m_final <= ncol(svd_result$u)) {
    B_reconstructor <- svd_result$u[, 1:m_final, drop = FALSE]
    
    # Compute coordinates (for compatibility with manifold output)
    if (m_final == 1) {
      Phi_coords <- svd_result$v[, 1, drop = FALSE] * svd_result$d[1]
    } else {
      Phi_coords <- svd_result$v[, 1:m_final, drop = FALSE] %*% 
                    diag(svd_result$d[1:m_final], m_final, m_final)
    }
  } else {
    # Fallback to constant
    warning("Using constant basis as ultimate fallback")
    B_reconstructor <- matrix(1/sqrt(p), p, 1)
    Phi_coords <- matrix(1, N, 1)
    m_final <- 1
  }
  
  return(list(
    B_reconstructor_matrix = B_reconstructor,
    Phi_coords_matrix = Phi_coords,
    eigenvalues_S_vector = c(1, var_explained),  # Add trivial eigenvalue
    m_final_dim = m_final,
    m_auto_selected_dim = m_auto,
    method_used = "PCA",
    variance_explained = cum_var[m_final]
  ))
}


#' Enhanced Manifold Basis Reconstructor with Fallback
#'
#' Robust version of get_manifold_basis_reconstructor_core with PCA fallback
#'
#' @param S_markov_matrix N x N Markov transition matrix
#' @param L_library_matrix p x N HRF library matrix
#' @param m_manifold_dim_target Target manifold dimensionality
#' @param m_manifold_dim_min_variance Minimum variance threshold
#' @param fallback_to_pca Whether to fall back to PCA on failure
#' @return List with reconstructor matrix and metadata
#' @export
get_manifold_basis_reconstructor_robust <- function(S_markov_matrix,
                                                   L_library_matrix,
                                                   m_manifold_dim_target,
                                                   m_manifold_dim_min_variance = 0.95,
                                                   fallback_to_pca = TRUE) {
  
  # First check library quality
  quality <- check_hrf_library_quality(L_library_matrix)
  
  if (!quality$is_good_quality) {
    warning("HRF library has quality issues:")
    if (!is.null(quality$has_duplicates) && !is.na(quality$has_duplicates) && quality$has_duplicates) {
      warning(sprintf("  - Found %d near-duplicate HRFs", quality$n_duplicates))
    }
    if (quality$is_ill_conditioned) {
      warning(sprintf("  - Library is ill-conditioned (condition number: %.2e)", 
                     quality$condition_number))
    }
    if (quality$has_degenerate_hrfs) {
      warning(sprintf("  - Found %d zero and %d constant HRFs", 
                     quality$n_zero_hrfs, quality$n_constant_hrfs))
    }
    if (quality$is_low_diversity) {
      warning("  - Low diversity in HRF shapes")
    }
  }
  
  # Try standard manifold construction
  result <- tryCatch({
    
    # Check if S_markov is degenerate
    if (is.matrix(S_markov_matrix)) {
      S_condition <- kappa(S_markov_matrix)
    } else {
      # For sparse matrices, check a subset
      S_sample <- as.matrix(S_markov_matrix[1:min(10, nrow(S_markov_matrix)), 
                                            1:min(10, ncol(S_markov_matrix))])
      S_condition <- kappa(S_sample)
    }
    
    if (S_condition > 1e12) {
      warning(sprintf("Markov matrix is poorly conditioned (kappa = %.2e)", S_condition))
      if (fallback_to_pca) {
        stop("Triggering PCA fallback due to poor conditioning")
      }
    }
    
    # Call original function
    result <- get_manifold_basis_reconstructor_core(
      S_markov_matrix = S_markov_matrix,
      L_library_matrix = L_library_matrix,
      m_manifold_dim_target = m_manifold_dim_target,
      m_manifold_dim_min_variance = m_manifold_dim_min_variance
    )
    
    # Check if result is reasonable
    if (any(is.na(result$B_reconstructor_matrix)) || 
        any(is.infinite(result$B_reconstructor_matrix))) {
      stop("Manifold construction produced invalid results")
    }
    
    # Add quality metrics
    result$library_quality <- quality
    result$method_used <- "diffusion_map"
    
    result
    
  }, error = function(e) {
    if (fallback_to_pca) {
      warning(sprintf("Manifold construction failed: %s", e$message))
      warning("Falling back to PCA-based approach")
      
      # Use PCA fallback
      pca_result <- compute_pca_fallback(
        L_library_matrix = L_library_matrix,
        m_target = m_manifold_dim_target,
        min_variance = m_manifold_dim_min_variance
      )
      
      pca_result$library_quality <- quality
      pca_result$original_error <- e$message
      
      return(pca_result)
    } else {
      stop(e)
    }
  })
  
  return(result)
}


#' Apply Physiological HRF Constraints
#'
#' Ensures HRF shapes are physiologically plausible
#'
#' @param hrf_matrix p x V matrix of HRF shapes (timepoints x voxels)
#' @param TR Time repetition in seconds
#' @param peak_range Acceptable peak time range in seconds (default c(2, 10))
#' @param enforce_positive Whether to enforce positive post-peak values
#' @param project_to_plausible Whether to project to plausible subspace
#' @return List with constrained HRFs and quality metrics
#' @export
apply_hrf_physiological_constraints <- function(hrf_matrix,
                                               TR = 2,
                                               peak_range = c(2, 10),
                                               enforce_positive = TRUE,
                                               project_to_plausible = TRUE) {
  
  p <- nrow(hrf_matrix)
  V <- ncol(hrf_matrix)
  time_vec <- (0:(p-1)) * TR
  
  # Initialize output
  hrf_constrained <- hrf_matrix
  quality_metrics <- matrix(NA, 4, V)
  rownames(quality_metrics) <- c("peak_time", "is_plausible", "adjustment_made", "integral")
  
  for (v in 1:V) {
    hrf <- hrf_matrix[, v]
    
    # Skip zero HRFs
    if (all(abs(hrf) < .Machine$double.eps)) {
      quality_metrics["is_plausible", v] <- FALSE
      next
    }
    
    # Find peak
    peak_idx <- which.max(hrf)
    peak_time <- time_vec[peak_idx]
    quality_metrics["peak_time", v] <- peak_time
    
    # Check if peak is in reasonable range
    peak_ok <- peak_time >= peak_range[1] && peak_time <= peak_range[2]
    
    # Check if integral is positive (net positive response)
    integral <- sum(hrf)
    quality_metrics["integral", v] <- integral
    integral_ok <- integral > 0
    
    # Check post-peak positivity (after initial undershoot)
    post_peak_ok <- TRUE
    if (peak_idx < p - 5) {  # Need at least 5 points after peak
      # Allow for initial undershoot but check late response
      late_response <- hrf[(peak_idx + 5):p]
      post_peak_ok <- mean(late_response) > -0.1 * max(hrf)
    }
    
    quality_metrics["is_plausible", v] <- peak_ok && integral_ok && post_peak_ok
    quality_metrics["adjustment_made", v] <- 0
    
    # Apply corrections if needed
    if (project_to_plausible && !quality_metrics["is_plausible", v]) {
      
      # Correct peak time if needed
      if (!peak_ok) {
        if (peak_time < peak_range[1]) {
          # Shift HRF to right
          shift_samples <- ceiling((peak_range[1] - peak_time) / TR)
          hrf_new <- c(rep(0, shift_samples), hrf[1:(p - shift_samples)])
        } else {
          # Shift HRF to left
          shift_samples <- ceiling((peak_time - peak_range[2]) / TR)
          hrf_new <- c(hrf[(shift_samples + 1):p], rep(0, shift_samples))
        }
        hrf_constrained[, v] <- hrf_new
        quality_metrics["adjustment_made", v] <- 1
      }
      
      # Ensure positive integral
      if (!integral_ok) {
        # Add small positive offset
        hrf_constrained[, v] <- hrf_constrained[, v] + 0.01
        quality_metrics["adjustment_made", v] <- 1
      }
      
      # Fix negative late response
      if (enforce_positive && !post_peak_ok && peak_idx < p - 5) {
        late_idx <- (peak_idx + 5):p
        hrf_constrained[late_idx, v] <- pmax(hrf_constrained[late_idx, v], 0)
        quality_metrics["adjustment_made", v] <- 1
      }
    }
  }
  
  # Compute overall reasonableness score
  n_plausible <- sum(quality_metrics["is_plausible", ], na.rm = TRUE)
  n_adjusted <- sum(quality_metrics["adjustment_made", ], na.rm = TRUE)
  
  return(list(
    hrf_constrained = hrf_constrained,
    quality_metrics = quality_metrics,
    percent_plausible = 100 * n_plausible / V,
    percent_adjusted = 100 * n_adjusted / V
  ))
}


#' Compute HRF Reasonableness Score
#'
#' Computes a score indicating how physiologically reasonable an HRF is
#'
#' @param hrf_vector Single HRF time course
#' @param TR Time repetition in seconds
#' @return Reasonableness score between 0 and 1
#' @keywords internal
compute_hrf_reasonableness <- function(hrf_vector, TR = 2) {
  
  p <- length(hrf_vector)
  time_vec <- (0:(p-1)) * TR
  
  # Zero HRF gets score 0
  if (all(abs(hrf_vector) < .Machine$double.eps)) {
    return(0)
  }
  
  # Normalize for comparison
  hrf_norm <- hrf_vector / max(abs(hrf_vector))
  
  # Score components
  scores <- numeric(5)
  
  # 1. Peak time score (best around 5s)
  peak_idx <- which.max(hrf_norm)
  peak_time <- time_vec[peak_idx]
  scores[1] <- exp(-0.5 * ((peak_time - 5) / 2)^2)  # Gaussian centered at 5s
  
  # 2. Peak width score (not too narrow or wide)
  half_max <- max(hrf_norm) / 2
  above_half <- which(hrf_norm > half_max)
  if (length(above_half) > 1) {
    fwhm <- (max(above_half) - min(above_half)) * TR
    scores[2] <- exp(-0.5 * ((fwhm - 6) / 3)^2)  # Ideal FWHM around 6s
  } else {
    scores[2] <- 0
  }
  
  # 3. Undershoot score (should have mild undershoot)
  if (peak_idx < p - 5) {
    undershoot_region <- hrf_norm[(peak_idx + 3):min(peak_idx + 10, p)]
    undershoot_depth <- -min(undershoot_region)
    # Ideal undershoot around 20% of peak
    scores[3] <- exp(-2 * abs(undershoot_depth - 0.2))
  } else {
    scores[3] <- 0.5  # Can't evaluate
  }
  
  # 4. Return to baseline score
  if (p > 15) {
    late_response <- mean(abs(hrf_norm[(p-5):p]))
    scores[4] <- exp(-10 * late_response)  # Should be near zero
  } else {
    scores[4] <- 0.5
  }
  
  # 5. Smoothness score (not too jagged)
  diff2 <- diff(diff(hrf_norm))
  roughness <- sum(diff2^2) / (p - 2)
  scores[5] <- exp(-5 * roughness)
  
  # Weighted average
  weights <- c(0.3, 0.2, 0.2, 0.15, 0.15)
  overall_score <- sum(scores * weights)
  
  return(overall_score)
}


