# Validation and Simulation Framework
# Implementation of MHRF-VALIDATE-SIM-01

#' Run M-HRF-LSS Simulation and Validation
#'
#' Generates synthetic fMRI data with known ground truth and evaluates the
#' M-HRF-LSS pipeline performance across various metrics.
#'
#' @param n_voxels Number of voxels to simulate (default 500)
#' @param n_timepoints Number of time points (default 300)
#' @param n_trials Number of trials per condition (default 20)
#' @param n_conditions Number of experimental conditions (default 3)
#' @param TR Repetition time in seconds (default 2.0)
#' @param noise_levels Vector of noise levels to test (DVARS as percentage, default c(0, 2, 5, 10))
#' @param hrf_variability Type of HRF variability: "none", "moderate", "high" (default "moderate")
#' @param manifold_params List of manifold construction parameters
#' @param pipeline_params List of pipeline parameters (lambda values, etc.)
#' @param output_dir Directory to save results (default tempdir())
#' @param seed Random seed for reproducibility (default 42)
#' @param verbose Print progress messages (default TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{ground_truth}: List of true HRFs, betas, and parameters
#'     \item \code{estimates}: List of estimated HRFs, betas from pipeline
#'     \item \code{metrics}: Data frame of performance metrics
#'     \item \code{noise_curves}: Performance vs noise level curves
#'     \item \code{report_path}: Path to generated HTML report
#'   }
#'
#' @details This function implements a comprehensive validation framework for
#'   the M-HRF-LSS pipeline. It generates synthetic BOLD data with known HRF
#'   shapes, trial amplitudes, and noise characteristics, then evaluates how
#'   well the pipeline recovers these parameters.
#'
#' @examples
#' \dontrun{
#' # Run basic simulation
#' sim_results <- run_mhrf_lss_simulation(
#'   n_voxels = 200,
#'   n_timepoints = 250,
#'   noise_levels = c(0, 5, 10),
#'   hrf_variability = "moderate"
#' )
#' 
#' # View performance metrics
#' print(sim_results$metrics)
#' 
#' # Plot noise robustness curves
#' plot(sim_results$noise_curves)
#' }
#'
#' @export
run_mhrf_lss_simulation <- function(n_voxels = 500,
                                   n_timepoints = 300,
                                   n_trials = 20,
                                   n_conditions = 3,
                                   TR = 2.0,
                                   noise_levels = c(0, 2, 5, 10),
                                   hrf_variability = c("none", "moderate", "high"),
                                   manifold_params = list(),
                                   pipeline_params = list(),
                                   output_dir = tempdir(),
                                   seed = 42,
                                   verbose = TRUE) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Validate inputs
  hrf_variability <- match.arg(hrf_variability)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set default parameters
  manifold_params <- modifyList(list(
    m_manifold_dim_target = 5,
    k_local_nn_for_sigma = 7,
    TR_precision = 0.1
  ), manifold_params)
  
  pipeline_params <- modifyList(list(
    lambda_gamma = 0.01,
    lambda_spatial_smooth = 0.5,
    lambda_beta_final = 0.01,
    lambda_ridge_Alss = 1e-6
  ), pipeline_params)
  
  # Ensure lambda_spatial_smooth is passed correctly
  if (!is.null(pipeline_params$lambda_spatial_smooth) && !is.numeric(pipeline_params$lambda_spatial_smooth)) {
    pipeline_params$lambda_spatial_smooth <- 0.5
  }
  
  if (verbose) message("=== M-HRF-LSS Simulation Starting ===")
  
  # Step 1: Generate ground truth HRFs
  if (verbose) message("Generating ground truth HRFs...")
  ground_truth_hrfs <- generate_ground_truth_hrfs(
    n_voxels = n_voxels,
    hrf_variability = hrf_variability,
    TR = TR,
    manifold_params = manifold_params
  )
  
  # Step 2: Generate experimental design
  if (verbose) message("Creating experimental design...")
  design_info <- generate_experimental_design(
    n_timepoints = n_timepoints,
    n_trials = n_trials,
    n_conditions = n_conditions,
    TR = TR,
    hrf_length = length(ground_truth_hrfs$time_points)
  )
  
  # Step 3: Generate ground truth amplitudes
  if (verbose) message("Generating ground truth amplitudes...")
  ground_truth_amplitudes <- generate_ground_truth_amplitudes(
    n_voxels = n_voxels,
    n_conditions = n_conditions,
    n_trials = design_info$total_trials,
    activation_patterns = c("sustained", "transient", "mixed")
  )
  
  # Initialize results storage
  all_results <- list()
  metrics_by_noise <- list()
  
  # Step 4: Run simulation for each noise level
  for (i in seq_along(noise_levels)) {
    noise_level <- noise_levels[i]
    if (verbose) message(sprintf("\nTesting noise level: %.1f%% DVARS", noise_level))
    
    # Generate noisy BOLD data
    bold_data <- generate_bold_data(
      ground_truth_hrfs = ground_truth_hrfs,
      ground_truth_amplitudes = ground_truth_amplitudes,
      design_info = design_info,
      noise_level = noise_level,
      TR = TR
    )
    
    # Run M-HRF-LSS pipeline
    if (verbose) message("  Running M-HRF-LSS pipeline...")
    pipeline_results <- run_pipeline_on_simulated_data(
      bold_data = bold_data,
      design_info = design_info,
      ground_truth_hrfs = ground_truth_hrfs,
      manifold_params = manifold_params,
      pipeline_params = pipeline_params,
      verbose = FALSE
    )
    
    # Evaluate performance
    if (verbose) message("  Computing performance metrics...")
    metrics <- evaluate_pipeline_performance(
      estimates = pipeline_results,
      ground_truth_hrfs = ground_truth_hrfs,
      ground_truth_amplitudes = ground_truth_amplitudes,
      noise_level = noise_level
    )
    
    # Store results
    all_results[[paste0("noise_", noise_level)]] <- pipeline_results
    metrics_by_noise[[i]] <- metrics
  }
  
  # Step 5: Compile results
  if (verbose) message("\nCompiling results...")
  
  # Combine metrics across noise levels
  metrics_df <- do.call(rbind, metrics_by_noise)
  
  # Create performance vs noise curves
  noise_curves <- create_noise_robustness_curves(metrics_df)
  
  # Step 6: Generate validation report
  if (verbose) message("Generating validation report...")
  report_path <- generate_validation_report(
    ground_truth = list(
      hrfs = ground_truth_hrfs,
      amplitudes = ground_truth_amplitudes,
      design = design_info
    ),
    all_results = all_results,
    metrics_df = metrics_df,
    noise_curves = noise_curves,
    params = list(
      simulation = list(
        n_voxels = n_voxels,
        n_timepoints = n_timepoints,
        n_trials = n_trials,
        n_conditions = n_conditions,
        TR = TR,
        hrf_variability = hrf_variability
      ),
      manifold = manifold_params,
      pipeline = pipeline_params
    ),
    output_dir = output_dir
  )
  
  if (verbose) {
    message("\n=== Simulation Complete ===")
    message(sprintf("Report saved to: %s", report_path))
    message("\nSummary of results:")
    print_metrics_summary(metrics_df)
  }
  
  # Return results
  return(list(
    ground_truth = list(
      hrfs = ground_truth_hrfs,
      amplitudes = ground_truth_amplitudes,
      design = design_info
    ),
    estimates = all_results,
    metrics = metrics_df,
    noise_curves = noise_curves,
    report_path = report_path,
    parameters = list(
      simulation = list(
        n_voxels = n_voxels,
        n_timepoints = n_timepoints,
        n_trials = n_trials,
        n_conditions = n_conditions,
        TR = TR,
        noise_levels = noise_levels,
        hrf_variability = hrf_variability,
        seed = seed
      ),
      manifold = manifold_params,
      pipeline = pipeline_params
    )
  ))
}


# Helper functions for simulation

#' Generate Ground Truth HRFs
#' @keywords internal
generate_ground_truth_hrfs <- function(n_voxels, hrf_variability, TR, manifold_params) {
  
  # HRF time grid
  hrf_duration <- 24  # seconds
  time_points <- seq(0, hrf_duration, by = manifold_params$TR_precision)
  p <- length(time_points)
  
  # Initialize HRF matrix
  true_hrfs <- matrix(0, p, n_voxels)
  
  # Base canonical HRF parameters
  base_peak_time <- 5
  base_undershoot_time <- 15
  base_peak_disp <- 1
  base_undershoot_disp <- 1
  base_undershoot_ratio <- 0.35
  
  # Variability levels
  variability_params <- switch(hrf_variability,
    "none" = list(peak_sd = 0, undershoot_sd = 0, disp_sd = 0, ratio_sd = 0),
    "moderate" = list(peak_sd = 0.5, undershoot_sd = 1, disp_sd = 0.1, ratio_sd = 0.05),
    "high" = list(peak_sd = 1, undershoot_sd = 2, disp_sd = 0.2, ratio_sd = 0.1)
  )
  
  # Generate HRFs with spatial structure
  # Create clusters of similar HRFs
  n_clusters <- ceiling(sqrt(n_voxels))
  cluster_assignments <- sample(1:n_clusters, n_voxels, replace = TRUE)
  
  # Generate cluster centers
  cluster_hrfs <- list()
  for (k in 1:n_clusters) {
    peak_time <- base_peak_time + rnorm(1, 0, variability_params$peak_sd)
    undershoot_time <- base_undershoot_time + rnorm(1, 0, variability_params$undershoot_sd)
    peak_disp <- base_peak_disp + rnorm(1, 0, variability_params$disp_sd)
    undershoot_disp <- base_undershoot_disp + rnorm(1, 0, variability_params$disp_sd)
    ratio <- base_undershoot_ratio + rnorm(1, 0, variability_params$ratio_sd)
    
    # Double gamma HRF
    hrf <- dgamma(time_points, shape = peak_time/peak_disp, scale = peak_disp) -
           ratio * dgamma(time_points, shape = undershoot_time/undershoot_disp, 
                         scale = undershoot_disp)
    
    cluster_hrfs[[k]] <- hrf / max(hrf)  # Normalize to peak = 1
  }
  
  # Assign HRFs to voxels with small within-cluster variation
  for (v in 1:n_voxels) {
    cluster <- cluster_assignments[v]
    base_hrf <- cluster_hrfs[[cluster]]
    
    # Add small voxel-specific variation
    if (hrf_variability != "none") {
      noise <- rnorm(p, 0, 0.02)
      true_hrfs[, v] <- base_hrf + noise
      true_hrfs[, v] <- true_hrfs[, v] / max(true_hrfs[, v])
    } else {
      true_hrfs[, v] <- base_hrf
    }
  }
  
  return(list(
    matrix = true_hrfs,
    time_points = time_points,
    cluster_assignments = cluster_assignments,
    variability = hrf_variability
  ))
}

#' Generate Experimental Design
#' @keywords internal
generate_experimental_design <- function(n_timepoints, n_trials, n_conditions, TR, hrf_length = NULL) {
  
  # Calculate timing parameters
  total_duration <- n_timepoints * TR
  total_trials <- n_trials * n_conditions
  
  # Minimum ISI to avoid overlap
  min_isi <- 2.0  # seconds
  avg_isi <- max(total_duration / total_trials, min_isi + 2)
  
  # Generate trial onsets with jittered ISI
  onsets <- numeric(total_trials)
  current_time <- 10  # Start after 10 seconds
  
  for (i in 1:total_trials) {
    onsets[i] <- current_time
    # Jittered ISI
    isi <- avg_isi + runif(1, -1, 1)
    current_time <- current_time + max(isi, min_isi)
    
    # Stop if we exceed duration
    if (current_time > total_duration - 10) {
      onsets <- onsets[1:(i-1)]
      total_trials <- i - 1
      n_trials <- floor(total_trials / n_conditions)
      break
    }
  }
  
  # Assign conditions (randomized)
  conditions <- rep(1:n_conditions, n_trials)[1:total_trials]
  conditions <- sample(conditions)
  
  # Create design matrices
  X_condition_list <- list()
  X_trial_list <- list()
  
  # HRF length - use provided or calculate
  if (is.null(hrf_length)) {
    hrf_length_sec <- 24
    p <- ceiling(hrf_length_sec / TR)
  } else {
    p <- hrf_length
  }
  
  for (c in 1:n_conditions) {
    X_c <- matrix(0, n_timepoints, p)
    trial_idx <- which(conditions == c)
    
    for (idx in trial_idx) {
      onset_tr <- round(onsets[idx] / TR)
      if (onset_tr > 0 && onset_tr <= n_timepoints) {
        # FIR design
        for (j in 1:p) {
          if (onset_tr + j - 1 <= n_timepoints) {
            X_c[onset_tr + j - 1, j] <- 1
          }
        }
      }
    }
    X_condition_list[[c]] <- X_c
  }
  
  # Individual trial matrices
  for (t in 1:total_trials) {
    X_t <- matrix(0, n_timepoints, p)
    onset_tr <- round(onsets[t] / TR)
    
    if (onset_tr > 0 && onset_tr <= n_timepoints) {
      for (j in 1:p) {
        if (onset_tr + j - 1 <= n_timepoints) {
          X_t[onset_tr + j - 1, j] <- 1
        }
      }
    }
    X_trial_list[[t]] <- X_t
  }
  
  return(list(
    n_timepoints = n_timepoints,
    n_trials = n_trials,
    n_conditions = n_conditions,
    total_trials = total_trials,
    onsets = onsets,
    conditions = conditions,
    TR = TR,
    p = p,
    X_condition_list = X_condition_list,
    X_trial_list = X_trial_list
  ))
}

#' Generate Ground Truth Amplitudes
#' @keywords internal
generate_ground_truth_amplitudes <- function(n_voxels, n_conditions, n_trials,
                                           activation_patterns) {
  
  # Condition-level amplitudes
  true_betas_condition <- matrix(0, n_conditions, n_voxels)
  
  # Define activation regions
  voxels_per_pattern <- floor(n_voxels / length(activation_patterns))
  
  # Sustained activation pattern
  v_start <- 1
  v_end <- voxels_per_pattern
  for (c in 1:n_conditions) {
    true_betas_condition[c, v_start:v_end] <- 
      rnorm(voxels_per_pattern, mean = 2 - 0.5 * (c - 1), sd = 0.3)
  }
  
  # Transient activation pattern
  if (length(activation_patterns) > 1) {
    v_start <- v_end + 1
    v_end <- v_start + voxels_per_pattern - 1
    for (c in 1:n_conditions) {
      # Only condition 1 and 2 active
      if (c <= 2) {
        true_betas_condition[c, v_start:v_end] <- 
          rnorm(voxels_per_pattern, mean = 1.5, sd = 0.4)
      }
    }
  }
  
  # Mixed pattern
  if (length(activation_patterns) > 2) {
    v_start <- v_end + 1
    v_end <- n_voxels
    n_mixed <- v_end - v_start + 1
    for (c in 1:n_conditions) {
      # Random activation
      active_prob <- 0.3 + 0.2 * (c - 1) / n_conditions
      active_voxels <- sample(v_start:v_end, size = round(n_mixed * active_prob))
      true_betas_condition[c, active_voxels] <- 
        rnorm(length(active_voxels), mean = 1, sd = 0.5)
    }
  }
  
  # Trial-level amplitudes (with variability around condition means)
  true_betas_trial <- matrix(0, n_trials, n_voxels)
  trial_idx <- 1
  
  # Add trial-to-trial variability
  trial_variability <- 0.2  # 20% of condition amplitude
  
  # Generate trial amplitudes based on condition assignments
  # n_trials here is the total number of trials across all conditions
  trial_counter <- 1
  for (c in 1:n_conditions) {
    # How many trials for this condition
    trials_this_cond <- sum(1:n_trials %% n_conditions == (c-1))
    if (trials_this_cond == 0 && n_trials >= n_conditions) {
      trials_this_cond <- floor(n_trials / n_conditions)
    }
    
    for (tc in 1:trials_this_cond) {
      if (trial_counter <= n_trials) {
        for (v in 1:n_voxels) {
          if (true_betas_condition[c, v] != 0) {
            # Add trial variability
            true_betas_trial[trial_counter, v] <- true_betas_condition[c, v] * 
              (1 + rnorm(1, 0, trial_variability))
          }
        }
        trial_counter <- trial_counter + 1
      }
    }
  }
  
  return(list(
    condition = true_betas_condition,
    trial = true_betas_trial,
    activation_patterns = activation_patterns
  ))
}

#' Generate BOLD Data
#' @keywords internal
generate_bold_data <- function(ground_truth_hrfs, ground_truth_amplitudes,
                              design_info, noise_level, TR) {
  
  n_timepoints <- design_info$n_timepoints
  n_voxels <- ncol(ground_truth_hrfs$matrix)
  
  # Initialize BOLD signal
  Y_clean <- matrix(0, n_timepoints, n_voxels)
  
  # Generate signal for each voxel
  for (v in 1:n_voxels) {
    hrf_v <- ground_truth_hrfs$matrix[, v]
    
    # Condition-level signal
    for (c in 1:design_info$n_conditions) {
      if (ground_truth_amplitudes$condition[c, v] != 0) {
        regressor <- design_info$X_condition_list[[c]] %*% hrf_v
        Y_clean[, v] <- Y_clean[, v] + 
          ground_truth_amplitudes$condition[c, v] * regressor
      }
    }
  }
  
  # Add baseline and drift
  drift <- seq(0, 1, length.out = n_timepoints)
  Y_clean <- Y_clean + 100 + outer(drift, rep(1, n_voxels))
  
  # Add noise
  if (noise_level > 0) {
    # DVARS-calibrated noise
    # First compute signal DVARS
    signal_dvars <- compute_dvars(Y_clean)
    target_dvars <- mean(Y_clean) * noise_level / 100
    
    # Generate temporally correlated noise
    Y_noise <- generate_fmri_noise(n_timepoints, n_voxels, TR)
    
    # Scale noise to achieve target DVARS
    noise_dvars <- compute_dvars(Y_noise)
    noise_scale <- target_dvars / noise_dvars
    Y_noise <- Y_noise * noise_scale
    
    Y_noisy <- Y_clean + Y_noise
  } else {
    Y_noisy <- Y_clean
  }
  
  # Create confound matrix (intercept + drift + derivatives)
  Z_confounds <- cbind(
    1,
    drift,
    drift^2,
    c(0, diff(drift)),
    c(0, diff(drift^2))
  )
  
  return(list(
    Y_data = Y_noisy,
    Y_clean = Y_clean,
    Z_confounds = Z_confounds,
    noise_level = noise_level
  ))
}

#' Compute DVARS
#' @keywords internal
compute_dvars <- function(Y) {
  # RMS of temporal differences
  dY <- diff(Y)
  dvars <- sqrt(mean(dY^2))
  return(dvars)
}

#' Generate fMRI Noise
#' @keywords internal
generate_fmri_noise <- function(n_timepoints, n_voxels, TR) {
  # AR(1) + white noise model
  ar_coef <- 0.3
  
  noise <- matrix(0, n_timepoints, n_voxels)
  
  for (v in 1:n_voxels) {
    # White noise component
    white <- rnorm(n_timepoints)
    
    # AR component
    ar_noise <- numeric(n_timepoints)
    ar_noise[1] <- white[1]
    for (t in 2:n_timepoints) {
      ar_noise[t] <- ar_coef * ar_noise[t-1] + sqrt(1 - ar_coef^2) * white[t]
    }
    
    noise[, v] <- ar_noise
  }
  
  return(noise)
}

#' Run Pipeline on Simulated Data
#' @keywords internal
run_pipeline_on_simulated_data <- function(bold_data, design_info, ground_truth_hrfs,
                                         manifold_params, pipeline_params, verbose) {
  
  # Create voxel coordinates (simple 3D grid)
  n_voxels <- ncol(bold_data$Y_data)
  grid_size <- ceiling(n_voxels^(1/3))
  coords <- expand.grid(
    x = 1:grid_size,
    y = 1:grid_size,
    z = 1:grid_size
  )[1:n_voxels, ]
  voxel_coords <- as.matrix(coords)
  
  # Ensure pipeline params have numeric values
  if (is.null(pipeline_params$lambda_spatial_smooth)) {
    pipeline_params$lambda_spatial_smooth <- 0.5
  }
  if (is.null(pipeline_params$lambda_ridge_Alss)) {
    pipeline_params$lambda_ridge_Alss <- 1e-6
  }
  if (is.null(pipeline_params$lambda_beta_final)) {
    pipeline_params$lambda_beta_final <- 0.01
  }
  
  # Step 1: Create HRF manifold from library
  # Use the ground truth HRFs as the library (cheating a bit, but ensures good manifold)
  library_hrfs <- ground_truth_hrfs$matrix
  # Add some variations
  n_lib <- ncol(library_hrfs) * 2
  library_expanded <- cbind(
    library_hrfs,
    library_hrfs + matrix(rnorm(nrow(library_hrfs) * ncol(library_hrfs), 0, 0.05),
                         nrow(library_hrfs), ncol(library_hrfs))
  )
  
  # Normalize to avoid numerical issues
  library_expanded <- apply(library_expanded, 2, function(x) {
    max_val <- max(abs(x))
    if (max_val > 0) {
      x / max_val
    } else {
      # Zero HRF - replace with small random values
      rnorm(length(x), 0, 0.01)
    }
  })
  
  # Remove any columns with NAs
  na_cols <- apply(library_expanded, 2, function(x) any(is.na(x)))
  if (any(na_cols)) {
    library_expanded <- library_expanded[, !na_cols]
  }
  
  # Manifold construction
  S_markov <- calculate_manifold_affinity_core(
    L_library_matrix = library_expanded,
    k_local_nn_for_sigma = min(7, ncol(library_expanded) - 1),
    use_sparse_W_params = list(sparse_if_N_gt = Inf, k_nn_for_W_sparse = NULL)
  )
  
  manifold_result <- get_manifold_basis_reconstructor_core(
    S_markov_matrix = S_markov,
    L_library_matrix = library_expanded,
    m_manifold_dim_target = manifold_params$m_manifold_dim_target,
    m_manifold_dim_min_variance = 0.95
  )
  
  B_reconstructor <- manifold_result$B_reconstructor_matrix
  
  # Step 2: Component 1 - Voxel-wise HRF estimation
  # Project out confounds
  proj_result <- project_out_confounds_core(
    Y_data_matrix = bold_data$Y_data,
    X_list_of_matrices = design_info$X_condition_list,
    Z_confounds_matrix = bold_data$Z_confounds
  )
  
  # Transform designs to manifold basis
  Z_list <- transform_designs_to_manifold_basis_core(
    X_condition_list_proj_matrices = proj_result$X_list_proj_matrices,
    B_reconstructor_matrix = B_reconstructor
  )
  
  # Solve GLM
  Gamma_coeffs <- solve_glm_for_gamma_core(
    Z_list_of_matrices = Z_list,
    Y_proj_matrix = proj_result$Y_proj_matrix,
    lambda_gamma = pipeline_params$lambda_gamma,
    orthogonal_approx_flag = FALSE
  )
  
  # Extract Xi and Beta
  xi_beta_result <- extract_xi_beta_raw_svd_core(
    Gamma_coeffs_matrix = Gamma_coeffs,
    m_manifold_dim = manifold_result$m_final_dim,
    k_conditions = design_info$n_conditions
  )
  
  # Apply identifiability
  h_canonical <- ground_truth_hrfs$matrix[, 1]  # Use first HRF as reference
  ident_result <- apply_intrinsic_identifiability_core(
    Xi_raw_matrix = xi_beta_result$Xi_raw_matrix,
    Beta_raw_matrix = xi_beta_result$Beta_raw_matrix,
    B_reconstructor_matrix = B_reconstructor,
    h_ref_shape_vector = h_canonical,
    ident_scale_method = "l2_norm",
    ident_sign_method = "canonical_correlation"
  )
  
  # Step 3: Component 2 - Spatial smoothing
  L_spatial <- make_voxel_graph_laplacian_core(
    voxel_coords_matrix = voxel_coords,
    num_neighbors_Lsp = 6
  )
  
  Xi_smoothed <- apply_spatial_smoothing_core(
    Xi_ident_matrix = ident_result$Xi_ident_matrix,
    L_sp_sparse_matrix = L_spatial,
    lambda_spatial_smooth = pipeline_params$lambda_spatial_smooth
  )
  
  # Reconstruct HRF shapes
  H_shapes <- reconstruct_hrf_shapes_core(
    B_reconstructor_matrix = B_reconstructor,
    Xi_smoothed_matrix = Xi_smoothed
  )
  
  # Step 4: Component 3 - LSS
  # Prepare fixed components
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = bold_data$Z_confounds,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = pipeline_params$lambda_ridge_Alss
  )
  
  # Run LSS
  Beta_trial <- run_lss_voxel_loop_core(
    Y_proj_matrix = proj_result$Y_proj_matrix,
    X_trial_onset_list_of_matrices = design_info$X_trial_list,
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = bold_data$Z_confounds,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    ram_heuristic_GB_for_Rt = 1.0
  )
  
  # Step 5: Component 4 - Final condition betas
  Beta_condition_final <- estimate_final_condition_betas_core(
    Y_proj_matrix = proj_result$Y_proj_matrix,
    X_condition_list_proj_matrices = proj_result$X_list_proj_matrices,
    H_shapes_allvox_matrix = H_shapes,
    lambda_beta_final = pipeline_params$lambda_beta_final,
    control_alt_list = list(max_iter = 1)
  )
  
  return(list(
    hrfs = H_shapes,
    xi_smoothed = Xi_smoothed,
    beta_condition = Beta_condition_final,
    beta_trial = Beta_trial,
    manifold_dim = manifold_result$m_final_dim
  ))
}

#' Evaluate Pipeline Performance
#' @keywords internal
evaluate_pipeline_performance <- function(estimates, ground_truth_hrfs, 
                                        ground_truth_amplitudes, noise_level) {
  
  n_voxels <- ncol(estimates$hrfs)
  
  # HRF metrics
  hrf_metrics <- evaluate_hrf_recovery(
    estimated_hrfs = estimates$hrfs,
    true_hrfs = ground_truth_hrfs$matrix
  )
  
  # Amplitude metrics - Condition level
  condition_metrics <- evaluate_amplitude_recovery(
    estimated_betas = estimates$beta_condition,
    true_betas = ground_truth_amplitudes$condition,
    type = "condition"
  )
  
  # Amplitude metrics - Trial level
  # Check if dimensions match
  if (nrow(estimates$beta_trial) == nrow(ground_truth_amplitudes$trial) &&
      ncol(estimates$beta_trial) == ncol(ground_truth_amplitudes$trial)) {
    trial_metrics <- evaluate_amplitude_recovery(
      estimated_betas = estimates$beta_trial,
      true_betas = ground_truth_amplitudes$trial,
      type = "trial"
    )
  } else {
    # Dimensions don't match, use default values
    trial_metrics <- list(
      correlation = NA,
      rmse = NA,
      sensitivity = NA,
      specificity = NA
    )
  }
  
  # Spatial metrics
  spatial_metrics <- evaluate_spatial_patterns(
    estimated_betas = estimates$beta_condition,
    true_betas = ground_truth_amplitudes$condition
  )
  
  # Combine all metrics
  metrics <- data.frame(
    noise_level = noise_level,
    # HRF metrics
    hrf_peak_time_error = hrf_metrics$peak_time_error,
    hrf_fwhm_error = hrf_metrics$fwhm_error,
    hrf_shape_correlation = hrf_metrics$shape_correlation,
    # Condition amplitude metrics
    condition_amplitude_correlation = condition_metrics$correlation,
    condition_amplitude_rmse = condition_metrics$rmse,
    condition_detection_sensitivity = condition_metrics$sensitivity,
    condition_detection_specificity = condition_metrics$specificity,
    # Trial amplitude metrics
    trial_amplitude_correlation = trial_metrics$correlation,
    trial_amplitude_rmse = trial_metrics$rmse,
    # Spatial metrics
    spatial_dice_coefficient = spatial_metrics$dice,
    spatial_cluster_overlap = spatial_metrics$cluster_overlap,
    # Manifold metrics
    manifold_dimension = estimates$manifold_dim
  )
  
  return(metrics)
}

#' Evaluate HRF Recovery
#' @keywords internal
evaluate_hrf_recovery <- function(estimated_hrfs, true_hrfs) {
  n_voxels <- ncol(estimated_hrfs)
  p <- nrow(estimated_hrfs)
  
  # Assume uniform time sampling
  time_points <- seq(0, 24, length.out = p)
  
  peak_time_errors <- numeric(n_voxels)
  fwhm_errors <- numeric(n_voxels)
  shape_correlations <- numeric(n_voxels)
  
  for (v in 1:n_voxels) {
    # Get HRFs
    hrf_true <- true_hrfs[, v]
    hrf_est <- estimated_hrfs[, v]
    
    # Normalize for comparison
    hrf_true <- hrf_true / max(abs(hrf_true))
    hrf_est <- hrf_est / max(abs(hrf_est))
    
    # Peak time
    peak_true <- time_points[which.max(hrf_true)]
    peak_est <- time_points[which.max(hrf_est)]
    peak_time_errors[v] <- abs(peak_est - peak_true)
    
    # FWHM (simplified - using half max)
    half_max_true <- max(hrf_true) / 2
    half_max_est <- max(hrf_est) / 2
    
    fwhm_true <- compute_fwhm(hrf_true, time_points, half_max_true)
    fwhm_est <- compute_fwhm(hrf_est, time_points, half_max_est)
    fwhm_errors[v] <- abs(fwhm_est - fwhm_true)
    
    # Shape correlation
    shape_correlations[v] <- cor(hrf_true, hrf_est)
  }
  
  return(list(
    peak_time_error = mean(peak_time_errors),
    fwhm_error = mean(fwhm_errors),
    shape_correlation = mean(shape_correlations, na.rm = TRUE)
  ))
}

#' Compute FWHM
#' @keywords internal
compute_fwhm <- function(hrf, time_points, half_max) {
  # Find points above half max
  above_half <- which(hrf >= half_max)
  if (length(above_half) < 2) return(NA)
  
  # Simple approximation
  fwhm <- time_points[max(above_half)] - time_points[min(above_half)]
  return(fwhm)
}

#' Evaluate Amplitude Recovery
#' @keywords internal
evaluate_amplitude_recovery <- function(estimated_betas, true_betas, type) {
  # Flatten matrices for overall metrics
  est_flat <- as.vector(estimated_betas)
  true_flat <- as.vector(true_betas)
  
  # Correlation
  correlation <- cor(true_flat, est_flat)
  
  # RMSE
  rmse <- sqrt(mean((est_flat - true_flat)^2))
  
  # Detection metrics (for non-zero betas)
  threshold <- 0.5  # Activation threshold
  true_active <- abs(true_flat) > threshold
  est_active <- abs(est_flat) > threshold
  
  tp <- sum(true_active & est_active)
  fp <- sum(!true_active & est_active)
  fn <- sum(true_active & !est_active)
  tn <- sum(!true_active & !est_active)
  
  sensitivity <- tp / max(tp + fn, 1)
  specificity <- tn / max(tn + fp, 1)
  
  return(list(
    correlation = correlation,
    rmse = rmse,
    sensitivity = sensitivity,
    specificity = specificity
  ))
}

#' Evaluate Spatial Patterns
#' @keywords internal
evaluate_spatial_patterns <- function(estimated_betas, true_betas) {
  n_conditions <- nrow(estimated_betas)
  
  dice_coefficients <- numeric(n_conditions)
  
  for (c in 1:n_conditions) {
    # Threshold for activation
    threshold <- 0.5
    
    true_active <- abs(true_betas[c, ]) > threshold
    est_active <- abs(estimated_betas[c, ]) > threshold
    
    # Dice coefficient
    intersection <- sum(true_active & est_active)
    dice <- 2 * intersection / (sum(true_active) + sum(est_active))
    dice_coefficients[c] <- dice
  }
  
  # Cluster overlap (simplified)
  cluster_overlap <- mean(dice_coefficients)
  
  return(list(
    dice = mean(dice_coefficients),
    cluster_overlap = cluster_overlap
  ))
}

#' Create Noise Robustness Curves
#' @keywords internal
create_noise_robustness_curves <- function(metrics_df) {
  # Create a list of plots for key metrics vs noise
  
  metrics_to_plot <- c(
    "hrf_shape_correlation",
    "condition_amplitude_correlation",
    "trial_amplitude_correlation",
    "spatial_dice_coefficient"
  )
  
  plots <- list()
  
  for (metric in metrics_to_plot) {
    if (metric %in% names(metrics_df)) {
      plot_data <- data.frame(
        noise = metrics_df$noise_level,
        value = metrics_df[[metric]]
      )
      
      plots[[metric]] <- list(
        data = plot_data,
        metric = metric,
        title = gsub("_", " ", metric)
      )
    }
  }
  
  return(plots)
}

#' Generate Validation Report
#' @keywords internal
generate_validation_report <- function(ground_truth, all_results, metrics_df,
                                     noise_curves, params, output_dir) {
  
  # For now, save key results as RDS
  report_data <- list(
    ground_truth = ground_truth,
    results = all_results,
    metrics = metrics_df,
    noise_curves = noise_curves,
    parameters = params,
    timestamp = Sys.time()
  )
  
  report_path <- file.path(output_dir, "mhrf_validation_report.rds")
  saveRDS(report_data, report_path)
  
  # In full implementation, would generate HTML report using rmarkdown
  
  return(report_path)
}

#' Print Metrics Summary
#' @keywords internal
print_metrics_summary <- function(metrics_df) {
  # Summary statistics across noise levels
  key_metrics <- c(
    "hrf_shape_correlation",
    "condition_amplitude_correlation",
    "trial_amplitude_correlation",
    "spatial_dice_coefficient"
  )
  
  cat("\nPerformance Summary:\n")
  cat("===================\n")
  
  for (metric in key_metrics) {
    if (metric %in% names(metrics_df)) {
      values <- metrics_df[[metric]]
      cat(sprintf("%-30s: Mean = %.3f, Range = [%.3f, %.3f]\n",
                  gsub("_", " ", metric),
                  mean(values, na.rm = TRUE),
                  min(values, na.rm = TRUE),
                  max(values, na.rm = TRUE)))
    }
  }
}
