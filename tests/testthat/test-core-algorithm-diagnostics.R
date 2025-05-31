# Diagnostic Tests for Core M-HRF-LSS Algorithm
# These tests verify mathematical correctness and algorithmic stability

test_that("M-HRF-LSS preserves signal reconstruction fidelity and manifold geometry", {
  # This test verifies that:
  # 1. The manifold embedding preserves HRF relationships
  # 2. Signal can be reconstructed accurately
  # 3. The method is stable under different SNR conditions
  
  set.seed(42)
  
  # Create ground truth HRF library with known structure
  p <- 30  # HRF length
  N <- 40  # Number of HRFs
  
  # Generate HRFs with systematic variations
  time_grid <- seq(0, 29, by = 1)
  
  # Create base HRF (double gamma)
  create_double_gamma <- function(peak_time, undershoot_ratio) {
    t <- time_grid
    # Positive gamma
    a1 <- peak_time
    b1 <- 1
    g1 <- (t/a1)^(a1/b1) * exp(-(t-a1)/b1)
    g1 <- g1 / max(g1)
    
    # Negative gamma (undershoot)
    a2 <- peak_time + 10
    b2 <- 1
    g2 <- (t/a2)^(a2/b2) * exp(-(t-a2)/b2)
    g2 <- g2 / max(g2) * undershoot_ratio
    
    hrf <- g1 - g2
    hrf / sum(abs(hrf))  # Normalize
  }
  
  # Generate library with systematic variations
  peak_times <- seq(4, 8, length.out = 8)
  undershoot_ratios <- seq(0, 0.5, length.out = 5)
  
  L_true <- matrix(0, p, N)
  idx <- 1
  for (pt in peak_times) {
    for (ur in undershoot_ratios) {
      if (idx <= N) {
        L_true[, idx] <- create_double_gamma(pt, ur)
        idx <- idx + 1
      }
    }
  }
  
  # TEST 1: Manifold preserves HRF geometry
  # Calculate true pairwise distances
  true_distances <- as.matrix(dist(t(L_true)))
  
  # Build manifold
  S_markov <- calculate_manifold_affinity_core(
    L_library_matrix = L_true,
    k_local_nn_for_sigma = 7
  )
  
  manifold <- get_manifold_basis_reconstructor_core(
    S_markov_matrix = S_markov,
    L_library_matrix = L_true,
    m_manifold_dim_target = 5,
    m_manifold_dim_min_variance = 0.95
  )
  
  # Check eigenvalue decay (should be smooth for good manifold)
  eigenvalues <- manifold$eigenvalues_S_vector[-1]  # Remove trivial eigenvalue
  eigenvalue_decay_rate <- diff(log(eigenvalues[1:5]))
  expect_true(all(eigenvalue_decay_rate < 0), 
              "Eigenvalues should decay monotonically")
  
  # Reconstruct HRFs from manifold coordinates
  # The reconstruction should be: L_approx = B * Phi'
  # But we need to get the manifold coordinates that correspond to our library
  # Since Phi_coords_matrix is N x m, we need: L_approx = B * Phi'
  # However, the correct reconstruction for the library itself is to project and reconstruct:
  
  # Use manifold coordinates for the library directly
  # Xi = transpose of Phi_coords_matrix (N x m -> m x N)
  Xi_library <- t(manifold$Phi_coords_matrix)

  # Now reconstruct HRFs
  L_reconstructed <- manifold$B_reconstructor_matrix %*% Xi_library
  
  # Check reconstruction error
  reconstruction_error <- norm(L_true - L_reconstructed, "F") / norm(L_true, "F")
  expect_lt(reconstruction_error, 0.1, 
            "HRF reconstruction error should be less than 10%")
  
  # Check that manifold preserves local neighborhoods
  # For each HRF, check if its k nearest neighbors are preserved
  k_check <- 5
  neighborhood_preservation <- numeric(N)
  
  for (i in 1:N) {
    # Original neighbors
    orig_neighbors <- order(true_distances[i, ])[2:(k_check+1)]
    
    # Manifold neighbors
    # Check if Phi_coords_matrix exists and has correct dimensions
    if (is.null(manifold$Phi_coords_matrix)) {
      # If not, compute from eigendecomposition
      manifold_coords <- manifold$B_reconstructor_matrix %*% diag(sqrt(manifold$eigenvalues_S_vector[2:(manifold$m_manifold_dim+1)]))
    } else {
      manifold_coords <- manifold$Phi_coords_matrix
    }
    
    # Ensure coordinates are N x m (HRFs as rows)
    if (nrow(manifold_coords) != N) {
      manifold_coords <- t(manifold_coords)
    }
    
    manifold_coords_dist <- as.matrix(dist(manifold_coords))
    manifold_neighbors <- order(manifold_coords_dist[i, ])[2:(k_check+1)]
    
    # Compute overlap
    neighborhood_preservation[i] <- length(intersect(orig_neighbors, manifold_neighbors)) / k_check
  }
  
  mean_preservation <- mean(neighborhood_preservation)
  expect_gt(mean_preservation, 0.7, 
            "Manifold should preserve at least 70% of local neighborhoods")
  
  # TEST 2: Signal reconstruction under different noise levels
  n_time <- 200
  n_voxels <- 50
  n_conditions <- 3
  n_trials_per_cond <- 10
  
  # Create realistic event design
  TR <- 2
  event_onsets <- sort(runif(n_conditions * n_trials_per_cond, 10, n_time * TR - 20))
  event_conditions <- rep(1:n_conditions, each = n_trials_per_cond)
  
  # Create design matrices
  X_conditions <- list()
  for (c in 1:n_conditions) {
    X_c <- matrix(0, n_time, p)
    trial_onsets <- event_onsets[event_conditions == c]
    
    for (onset in trial_onsets) {
      time_idx <- floor(onset / TR)
      if (time_idx + p <= n_time) {
        X_c[time_idx:(time_idx + p - 1), ] <- X_c[time_idx:(time_idx + p - 1), ] + diag(p)
      }
    }
    X_conditions[[c]] <- X_c
  }
  
  # Generate ground truth: Each voxel has a random HRF from our library
  true_hrf_indices <- sample(1:N, n_voxels, replace = TRUE)
  H_true <- L_true[, true_hrf_indices]
  
  # Generate true betas with some structure
  Beta_true <- matrix(0, n_conditions, n_voxels)
  for (v in 1:n_voxels) {
    # Create correlated betas across conditions
    base_activation <- runif(1, 0.5, 2)
    Beta_true[, v] <- base_activation * (1 + rnorm(n_conditions, 0, 0.3))
  }
  Beta_true[Beta_true < 0] <- 0
  
  # Generate signal
  Y_clean <- matrix(0, n_time, n_voxels)
  for (c in 1:n_conditions) {
    # X_conditions[[c]] is n_time x p
    # H_true is p x n_voxels  
    # Beta_true[c, ] is a vector of length n_voxels
    # We want: Y = X * H * diag(beta) where result is n_time x n_voxels
    Y_clean <- Y_clean + X_conditions[[c]] %*% (H_true * matrix(Beta_true[c, ], p, n_voxels, byrow = TRUE))
  }
  
  # Test multiple SNR levels
  snr_levels <- c(Inf, 2, 1, 0.5)  # Inf = no noise
  recovery_errors <- list()
  
  for (snr in snr_levels) {
    # Add noise
    if (is.finite(snr)) {
      signal_power <- mean(Y_clean^2)
      noise_power <- signal_power / snr
      Y_noisy <- Y_clean + matrix(rnorm(n_time * n_voxels, sd = sqrt(noise_power)), n_time, n_voxels)
    } else {
      Y_noisy <- Y_clean
    }
    
    # Add confounds
    Z_confounds <- cbind(1, poly(1:n_time, degree = 2))
    
    # Run M-HRF-LSS pipeline
    proj_result <- project_out_confounds_core(
      Y_data_matrix = Y_noisy,
      X_list_of_matrices = X_conditions,
      Z_confounds_matrix = Z_confounds
    )
    
    # Transform to manifold basis
    Z_list <- transform_designs_to_manifold_basis_core(
      X_condition_list_proj_matrices = proj_result$X_list_proj_matrices,
      B_reconstructor_matrix = manifold$B_reconstructor_matrix
    )
    
    # Solve for gamma
    Gamma_est <- solve_glm_for_gamma_core(
      Z_list_of_matrices = Z_list,
      Y_proj_matrix = proj_result$Y_proj_matrix,
      lambda_gamma = 0.01
    )
    
    # Extract Xi and Beta
    xi_beta <- extract_xi_beta_raw_svd_robust(
      Gamma_coeffs_matrix = Gamma_est,
      m_manifold_dim = manifold$m_manifold_dim,
      k_conditions = n_conditions
    )
    
    # Apply identifiability
    h_ref <- L_true[, 1]
    ident_result <- apply_intrinsic_identifiability_core(
      Xi_raw_matrix = xi_beta$Xi_raw_matrix,
      Beta_raw_matrix = xi_beta$Beta_raw_matrix,
      B_reconstructor_matrix = manifold$B_reconstructor_matrix,
      h_ref_shape_vector = h_ref,
      ident_scale_method = "l2_norm",
      ident_sign_method = "canonical_correlation"
    )
    
    # Reconstruct HRFs
    H_est <- reconstruct_hrf_shapes_core(
      B_reconstructor_matrix = manifold$B_reconstructor_matrix,
      Xi_smoothed_matrix = ident_result$Xi_ident_matrix
    )
    
    # Compute recovery errors
    # 1. HRF shape recovery (up to scale/sign)
    hrf_correlations <- numeric(n_voxels)
    for (v in 1:n_voxels) {
      if (sum(H_est[, v]^2) > 0 && sum(H_true[, v]^2) > 0) {
        hrf_correlations[v] <- abs(cor(H_est[, v], H_true[, v]))
      } else {
        hrf_correlations[v] <- 0
      }
    }
    
    # 2. Beta recovery (relative error)
    beta_error <- norm(ident_result$Beta_ident_matrix - Beta_true, "F") / norm(Beta_true, "F")
    
    recovery_errors[[as.character(snr)]] <- list(
      hrf_correlation = mean(hrf_correlations, na.rm = TRUE),
      beta_relative_error = beta_error,
      n_degenerate_voxels = sum(hrf_correlations == 0)
    )
  }
  
  # Verify recovery quality
  # With no noise, should have excellent recovery
  expect_gt(recovery_errors[["Inf"]]$hrf_correlation, 0.95,
            "HRF recovery should be >95% correlation with no noise")
  expect_lt(recovery_errors[["Inf"]]$beta_relative_error, 0.1,
            "Beta recovery error should be <10% with no noise")
  
  # With SNR=2, should still have good recovery
  expect_gt(recovery_errors[["2"]]$hrf_correlation, 0.8,
            "HRF recovery should be >80% correlation at SNR=2")
  expect_lt(recovery_errors[["2"]]$beta_relative_error, 0.3,
            "Beta recovery error should be <30% at SNR=2")
  
  # Verify graceful degradation
  snr_values <- c(Inf, 2, 1, 0.5)
  hrf_corrs <- sapply(recovery_errors, function(x) x$hrf_correlation)
  expect_true(all(diff(hrf_corrs) <= 0),
              "HRF recovery should degrade monotonically with decreasing SNR")
  
  # No voxels should be completely degenerate at reasonable SNR
  expect_equal(recovery_errors[["2"]]$n_degenerate_voxels, 0,
               "No voxels should be degenerate at SNR=2")
})


test_that("M-HRF-LSS trial-wise estimation is unbiased and efficient compared to standard LSS", {
  # This test verifies that:
  # 1. Trial-wise beta estimates are unbiased
  # 2. The Woodbury implementation is numerically equivalent to naive LSS
  # 3. Efficiency gains are achieved without sacrificing accuracy
  # 4. The method handles overlapping trials correctly
  
  set.seed(123)
  
  # Simulation parameters
  n_time <- 300
  n_voxels <- 20  # Smaller for detailed comparison
  n_conditions <- 2
  n_trials_per_cond <- 15
  TR <- 2
  p_hrf <- 25
  
  # Create true HRF (canonical double-gamma)
  t <- seq(0, (p_hrf-1)*TR, by = TR)
  hrf_true <- dgamma(t, shape = 6, scale = 1) - 0.35 * dgamma(t, shape = 16, scale = 1)
  hrf_true <- hrf_true / sum(abs(hrf_true))
  
  # Generate trial onsets with some close spacing to test overlap handling
  min_iti <- 4  # seconds
  max_iti <- 12
  
  trial_onsets <- numeric(n_conditions * n_trials_per_cond)
  current_time <- 10
  
  for (i in 1:(n_conditions * n_trials_per_cond)) {
    trial_onsets[i] <- current_time
    current_time <- current_time + runif(1, min_iti, max_iti)
  }
  
  trial_conditions <- rep(1:n_conditions, each = n_trials_per_cond)
  
  # Generate true trial-wise amplitudes with structured variance
  trial_betas_true <- matrix(0, n_conditions * n_trials_per_cond, n_voxels)
  
  for (v in 1:n_voxels) {
    # Each voxel has different condition effects
    cond_means <- runif(n_conditions, 0.5, 2)
    
    # Add trial-wise variability
    for (trial in 1:length(trial_onsets)) {
      cond <- trial_conditions[trial]
      trial_betas_true[trial, v] <- cond_means[cond] * exp(rnorm(1, 0, 0.3))
    }
  }
  
  # Create trial-wise design matrices
  X_trials <- list()
  for (trial in 1:length(trial_onsets)) {
    X_trial <- matrix(0, n_time, p_hrf)
    onset_idx <- round(trial_onsets[trial] / TR) + 1
    
    if (onset_idx + p_hrf - 1 <= n_time) {
      X_trial[onset_idx:(onset_idx + p_hrf - 1), ] <- diag(p_hrf)
    }
    X_trials[[trial]] <- X_trial
  }
  
  # Generate clean signal
  Y_clean <- matrix(0, n_time, n_voxels)
  for (trial in 1:length(trial_onsets)) {
    for (v in 1:n_voxels) {
      Y_clean[, v] <- Y_clean[, v] + X_trials[[trial]] %*% hrf_true * trial_betas_true[trial, v]
    }
  }
  
  # Add realistic confounds
  Z_confounds <- cbind(
    1,  # Intercept
    poly(1:n_time, degree = 2),  # Polynomial drift
    sin(2*pi*(1:n_time)/(n_time/5)),  # Physiological noise
    cos(2*pi*(1:n_time)/(n_time/5))
  )
  
  # Add measurement noise
  snr <- 1.5
  signal_power <- mean(apply(Y_clean, 2, var))
  noise_power <- signal_power / snr
  Y_noisy <- Y_clean + matrix(rnorm(n_time * n_voxels, sd = sqrt(noise_power)), n_time, n_voxels)
  
  # Add confound signal
  confound_weights <- matrix(rnorm(ncol(Z_confounds) * n_voxels), ncol(Z_confounds), n_voxels)
  Y_observed <- Y_noisy + Z_confounds %*% confound_weights
  
  # TEST 1: Woodbury LSS vs Naive LSS equivalence
  # Project out confounds from Y only (not from X matrices)
  # The Woodbury implementation expects unprojected X matrices
  P_confound <- diag(n_time) - Z_confounds %*% solve(crossprod(Z_confounds) + 1e-6 * diag(ncol(Z_confounds))) %*% t(Z_confounds)
  Y_proj <- P_confound %*% Y_observed
  
  # Prepare LSS components
  A_fixed <- Z_confounds
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = A_fixed,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = 1e-6
  )
  
  # For a subset of voxels, compare Woodbury vs naive implementation
  test_voxels <- sample(1:n_voxels, min(5, n_voxels))
  
  for (v in test_voxels) {
    # Get HRF for this voxel (using true HRF for this test)
    h_v <- hrf_true
    
    # Woodbury implementation: projected Y, unprojected X
    beta_woodbury <- run_lss_for_voxel_core(
      Y_proj_voxel_vector = Y_proj[, v],
      X_trial_onset_list_of_matrices = X_trials,  # UNPROJECTED
      H_shape_voxel_vector = h_v,
      A_lss_fixed_matrix = A_fixed,
      P_lss_matrix = lss_prep$P_lss_matrix,
      p_lss_vector = lss_prep$p_lss_vector
    )
    
    # Naive implementation (for comparison)
    beta_naive <- numeric(length(X_trials))
    for (trial in 1:length(X_trials)) {
      # Build full design matrix for this trial
      X_other_trials <- matrix(0, n_time, 0)
      for (j in setdiff(1:length(X_trials), trial)) {
        X_other_trials <- cbind(X_other_trials, X_trials[[j]] %*% h_v)
      }
      
      # Full design includes confounds
      X_full <- cbind(
        X_trials[[trial]] %*% h_v,  # Trial of interest
        X_other_trials,              # Other trials
        A_fixed                      # Confounds
      )
      
      # Solve via normal equations using UNPROJECTED Y
      XtX <- crossprod(X_full)
      XtX <- XtX + diag(1e-6, nrow(XtX))  # Ridge
      Xty <- crossprod(X_full, Y_observed[, v])  # Use original Y
      
      beta_all <- solve(XtX, Xty)
      beta_naive[trial] <- beta_all[1]  # First coefficient is trial of interest
    }
    
    # Compare results
    max_diff <- max(abs(beta_woodbury - beta_naive))
    expect_lt(max_diff, 1e-10,
              sprintf("Woodbury and naive LSS should agree to numerical precision (voxel %d)", v))
  }
  
  # TEST 2: Recovery of trial-wise betas
  # Run full M-HRF-LSS pipeline
  
  # For this test, use known HRF
  H_shapes <- matrix(rep(hrf_true, n_voxels), p_hrf, n_voxels)
  
  # Run LSS with correct data flow:
  # - Y_proj_matrix should be projected Y
  # - X_trial_onset_list should be UNPROJECTED
  beta_lss <- run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,  # Already computed above
    X_trial_onset_list_of_matrices = X_trials,  # UNPROJECTED
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = Z_confounds,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector
  )
  
  # Compute recovery metrics
  # 1. Bias: mean error
  bias <- mean(beta_lss - trial_betas_true)
  expect_lt(abs(bias), 0.05,
            "Trial-wise estimates should be unbiased (bias < 0.05)")
  
  # 2. Correlation with truth
  correlation <- cor(as.vector(beta_lss), as.vector(trial_betas_true))
  expect_gt(correlation, 0.85,
            "Trial-wise estimates should correlate >0.85 with truth at SNR=1.5")
  
  # 3. Variance explained
  ss_total <- sum((trial_betas_true - mean(trial_betas_true))^2)
  ss_residual <- sum((beta_lss - trial_betas_true)^2)
  r_squared <- 1 - ss_residual / ss_total
  expect_gt(r_squared, 0.7,
            "LSS should explain >70% variance in trial-wise betas")
  
  # TEST 3: Condition-wise aggregation preserves structure
  # Average trial betas within condition
  beta_condition_true <- matrix(0, n_conditions, n_voxels)
  beta_condition_lss <- matrix(0, n_conditions, n_voxels)
  
  for (c in 1:n_conditions) {
    trials_c <- which(trial_conditions == c)
    beta_condition_true[c, ] <- colMeans(trial_betas_true[trials_c, , drop = FALSE])
    beta_condition_lss[c, ] <- colMeans(beta_lss[trials_c, , drop = FALSE])
  }
  
  # Condition-level estimates should be even more accurate
  cond_correlation <- cor(as.vector(beta_condition_lss), as.vector(beta_condition_true))
  expect_gt(cond_correlation, 0.95,
            "Condition-level averages should correlate >0.95 with truth")
  
  # TEST 4: Overlapping trials are handled correctly
  # Find trials that overlap in time
  overlapping_pairs <- list()
  for (i in 1:(length(trial_onsets)-1)) {
    for (j in (i+1):length(trial_onsets)) {
      if (abs(trial_onsets[i] - trial_onsets[j]) < p_hrf * TR) {
        overlapping_pairs <- append(overlapping_pairs, list(c(i, j)))
      }
    }
  }
  
  expect_gt(length(overlapping_pairs), 0,
            "Test should include overlapping trials")
  
  # Check that overlapping trials can still be distinguished
  if (length(overlapping_pairs) > 0) {
    # Pick first overlapping pair
    pair <- overlapping_pairs[[1]]
    t1 <- pair[1]
    t2 <- pair[2]
    
    # Find voxels where these trials had different amplitudes
    amplitude_diff <- abs(trial_betas_true[t1, ] - trial_betas_true[t2, ])
    discriminative_voxels <- which(amplitude_diff > 0.5)
    
    if (length(discriminative_voxels) > 0) {
      v <- discriminative_voxels[1]
      
      # Check that LSS recovered different amplitudes
      lss_diff <- abs(beta_lss[t1, v] - beta_lss[t2, v])
      true_diff <- abs(trial_betas_true[t1, v] - trial_betas_true[t2, v])
      
      recovery_ratio <- lss_diff / true_diff
      expect_gt(recovery_ratio, 0.5,
                "LSS should recover at least 50% of amplitude difference for overlapping trials")
    }
  }
  
  # TEST 5: Efficiency check
  # Verify that trial-wise variance is reasonable
  for (v in 1:n_voxels) {
    # Within-condition variance
    for (c in 1:n_conditions) {
      trials_c <- which(trial_conditions == c)
      
      if (length(trials_c) > 1) {
        true_var <- var(trial_betas_true[trials_c, v])
        lss_var <- var(beta_lss[trials_c, v])
        
        # LSS variance should not be inflated by more than 3x (realistic for overlapping trials)
        if (true_var > 0.01) {  # Only check voxels with meaningful variance
          variance_inflation <- lss_var / true_var
          expect_lt(variance_inflation, 3,
                    sprintf("LSS variance inflation should be <3x (voxel %d, condition %d)", v, c))
        }
      }
    }
  }
})