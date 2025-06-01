# Benchmark Performance Script for M-HRF-LSS
# Compare performance across different fmrireg benchmark datasets

library(manifoldhrf)
library(fmrireg)
library(ggplot2)

# Function to evaluate performance on a benchmark dataset
evaluate_benchmark <- function(dataset_name, preset = "balanced", verbose = TRUE) {
  
  if (verbose) cat("\n=== Evaluating", dataset_name, "===\n")
  
  # Load data
  bm_data <- fmrireg:::load_benchmark_dataset(dataset_name)
  summary_info <- fmrireg:::get_benchmark_summary(dataset_name)
  
  Y_data <- bm_data$data
  event_list <- bm_data$event_list
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  
  # Timing
  start_time <- Sys.time()
  
  # Create design matrices
  p <- 25
  conditions <- unique(event_list$condition)
  k <- length(conditions)
  
  X_condition_list <- list()
  X_trial_list <- list()
  
  for (i in 1:k) {
    cond_events <- event_list[event_list$condition == conditions[i], ]
    X_cond <- matrix(0, n, p)
    
    for (j in 1:nrow(cond_events)) {
      onset_idx <- round(cond_events$onset[j] / 2) + 1
      duration_idx <- max(1, round(cond_events$duration[j] / 2))
      
      # Create trial matrix
      X_trial <- matrix(0, n, p)
      
      for (t in 0:(duration_idx - 1)) {
        if (onset_idx + t <= n - p + 1) {
          X_cond[(onset_idx + t):(onset_idx + t + p - 1), ] <- 
            X_cond[(onset_idx + t):(onset_idx + t + p - 1), ] + diag(p)
          
          if (i == 1 && j <= 10) {  # Store first 10 trials of first condition
            X_trial[(onset_idx + t):(onset_idx + t + p - 1), ] <- diag(p)
          }
        }
      }
      
      if (i == 1 && j <= 10) {
        X_trial_list[[length(X_trial_list) + 1]] <- X_trial
      }
    }
    X_condition_list[[i]] <- X_cond
  }
  
  # Get parameters
  params <- get_preset_params(preset, n_voxels = V)
  
  # Create HRF library
  L_library <- create_gamma_grid_library(
    n_samples = p,
    N_hrfs = 50,
    TR = 2,
    peak_range = c(3, 9),
    width_range = c(3, 8),
    normalize = TRUE
  )
  
  # Run M-HRF-LSS pipeline
  results <- list()
  
  tryCatch({
    # Component 0: Manifold construction
    manifold <- create_hrf_manifold(
      hrf_matrix = L_library,
      m_target = params$m_manifold_dim_target,
      k_nn = params$k_local_nn_for_sigma
    )
    
    # Component 1: Voxel-wise fit
    if (verbose) cat("Running voxel-wise HRF estimation...\n")
    
    # Project out confounds (none in this case)
    Y_clean <- Y_data
    X_clean <- X_condition_list
    
    # Transform to manifold basis
    XB_list <- transform_designs_to_manifold_basis_core(
      X_condition_list = X_clean,
      B_manifold_matrix = manifold$B_reconstructor
    )
    
    # Solve for gamma
    Gamma <- solve_glm_for_gamma_core(
      Y_data_matrix = Y_clean,
      XB_condition_list = XB_list,
      lambda_ridge = params$lambda_gamma
    )
    
    # Extract Xi and Beta
    svd_result <- extract_xi_beta_raw_svd_core(
      Gamma_coeffs_matrix = Gamma,
      m_manifold_dim = manifold$m_selected,
      k_conditions = k
    )
    
    # Apply identifiability
    ident_result <- apply_intrinsic_identifiability_core(
      Xi_raw_matrix = svd_result$Xi_raw_matrix,
      Beta_raw_matrix = svd_result$Beta_raw_matrix,
      Gamma_raw_matrix = Gamma
    )
    
    # Component 2: Spatial smoothing (minimal for benchmark)
    if (verbose) cat("Applying spatial smoothing...\n")
    Xi_smooth <- apply_spatial_smoothing_core(
      Xi_ident_matrix = ident_result$Xi_ident_matrix,
      L_sp_sparse_matrix = diag(V),  # No spatial structure in benchmark
      lambda_spatial_smooth = 0.01
    )
    
    # Get HRF shapes
    hrf_shapes <- reconstruct_hrf_shapes_core(
      Xi_smooth_matrix = Xi_smooth,
      B_reconstructor_matrix = manifold$B_reconstructor
    )
    
    # Component 3: Trial-wise LSS (sample)
    if (length(X_trial_list) > 0 && verbose) {
      cat("Running trial-wise LSS (sample)...\n")
      
      # Test on first voxel
      lss_result <- run_lss_for_voxel_corrected(
        y_voxel = Y_data[, 1],
        X_trial_list = X_trial_list,
        h_voxel = hrf_shapes[, 1],
        TR = 2,
        lambda = params$lambda_ridge_Alss
      )
      
      results$trial_betas_sample <- lss_result$beta_trials
    }
    
    results$success <- TRUE
    results$Xi_smooth <- Xi_smooth
    results$Beta_condition <- ident_result$Beta_ident_matrix
    results$hrf_shapes <- hrf_shapes
    results$manifold_method <- manifold$method_used
    
  }, error = function(e) {
    results$success <- FALSE
    results$error <- e$message
    if (verbose) cat("Error:", e$message, "\n")
  })
  
  end_time <- Sys.time()
  results$runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Evaluate performance metrics
  if (results$success) {
    metrics <- list()
    
    # 1. HRF recovery (if canonical)
    if (grepl("Canonical", dataset_name)) {
      hrf_canonical <- fmrireg::HRF_SPMG1(seq(0, by = 2, length.out = p))
      hrf_canonical <- hrf_canonical / sum(abs(hrf_canonical))
      
      hrf_correlations <- cor(hrf_canonical, results$hrf_shapes)
      metrics$hrf_recovery <- list(
        median_correlation = median(hrf_correlations),
        mean_correlation = mean(hrf_correlations),
        min_correlation = min(hrf_correlations)
      )
    }
    
    # 2. Reconstruction error
    Y_reconstructed <- matrix(0, n, V)
    for (i in 1:k) {
      X_conv <- matrix(0, n, V)
      for (v in 1:V) {
        # Convolve design with estimated HRF
        design_pad <- rbind(X_condition_list[[i]], matrix(0, p-1, p))
        hrf_v <- results$hrf_shapes[, v]
        
        for (t in 1:n) {
          X_conv[t, v] <- sum(design_pad[t:(t+p-1), ] %*% hrf_v)
        }
      }
      Y_reconstructed <- Y_reconstructed + X_conv * 
        matrix(results$Beta_condition[i, ], n, V, byrow = TRUE)
    }
    
    residuals <- Y_data - Y_reconstructed
    metrics$reconstruction <- list(
      rmse = sqrt(mean(residuals^2)),
      r_squared = 1 - sum(residuals^2) / sum((Y_data - mean(Y_data))^2),
      mean_voxel_r2 = mean(apply(1:V, 1, function(v) {
        1 - sum(residuals[, v]^2) / sum((Y_data[, v] - mean(Y_data[, v]))^2)
      }))
    )
    
    # 3. Condition separation (if multiple conditions)
    if (k > 1) {
      # F-statistic for condition differences
      beta_means <- rowMeans(results$Beta_condition)
      beta_vars <- apply(results$Beta_condition, 1, var)
      between_var <- var(beta_means) * V
      within_var <- mean(beta_vars)
      
      metrics$condition_separation <- list(
        f_statistic = between_var / within_var,
        effect_size = sqrt(between_var / (between_var + within_var))
      )
    }
    
    # 4. Computational efficiency
    metrics$efficiency <- list(
      voxels_per_second = V / results$runtime,
      total_time = results$runtime
    )
    
    results$metrics <- metrics
  }
  
  # Summary
  if (verbose) {
    cat("\nResults for", dataset_name, ":\n")
    cat("  Success:", results$success, "\n")
    if (results$success) {
      cat("  Runtime:", round(results$runtime, 2), "seconds\n")
      cat("  Manifold method:", results$manifold_method, "\n")
      if (!is.null(results$metrics$hrf_recovery)) {
        cat("  HRF recovery (median cor):", 
            round(results$metrics$hrf_recovery$median_correlation, 3), "\n")
      }
      cat("  Reconstruction R²:", 
          round(results$metrics$reconstruction$r_squared, 3), "\n")
      cat("  Mean voxel R²:", 
          round(results$metrics$reconstruction$mean_voxel_r2, 3), "\n")
    }
  }
  
  return(results)
}


# Run evaluation on all datasets
run_all_benchmarks <- function(presets = c("conservative", "balanced", "robust")) {
  
  datasets <- c(
    "BM_Canonical_HighSNR",
    "BM_Canonical_LowSNR", 
    "BM_HRF_Variability_AcrossVoxels",
    "BM_Trial_Amplitude_Variability",
    "BM_Complex_Realistic"
  )
  
  all_results <- list()
  
  for (preset in presets) {
    cat("\n\n========== PRESET:", preset, "==========\n")
    
    preset_results <- list()
    
    for (dataset in datasets) {
      result <- evaluate_benchmark(dataset, preset = preset)
      preset_results[[dataset]] <- result
    }
    
    all_results[[preset]] <- preset_results
  }
  
  return(all_results)
}


# Create summary plots
plot_benchmark_summary <- function(all_results) {
  
  # Extract metrics for plotting
  plot_data <- data.frame()
  
  for (preset in names(all_results)) {
    for (dataset in names(all_results[[preset]])) {
      result <- all_results[[preset]][[dataset]]
      
      if (result$success && !is.null(result$metrics)) {
        row <- data.frame(
          preset = preset,
          dataset = dataset,
          runtime = result$runtime,
          r_squared = result$metrics$reconstruction$r_squared,
          mean_voxel_r2 = result$metrics$reconstruction$mean_voxel_r2,
          stringsAsFactors = FALSE
        )
        
        if (!is.null(result$metrics$hrf_recovery)) {
          row$hrf_correlation <- result$metrics$hrf_recovery$median_correlation
        } else {
          row$hrf_correlation <- NA
        }
        
        plot_data <- rbind(plot_data, row)
      }
    }
  }
  
  # Create plots
  plots <- list()
  
  # 1. Reconstruction R² by dataset and preset
  plots$r2 <- ggplot(plot_data, aes(x = dataset, y = r_squared, fill = preset)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Reconstruction R² by Dataset and Preset",
         y = "R²", x = "Dataset") +
    ylim(0, 1)
  
  # 2. Runtime comparison
  plots$runtime <- ggplot(plot_data, aes(x = dataset, y = runtime, fill = preset)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Runtime by Dataset and Preset",
         y = "Runtime (seconds)", x = "Dataset")
  
  # 3. HRF recovery (for canonical datasets)
  hrf_data <- plot_data[!is.na(plot_data$hrf_correlation), ]
  if (nrow(hrf_data) > 0) {
    plots$hrf <- ggplot(hrf_data, aes(x = dataset, y = hrf_correlation, fill = preset)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "HRF Recovery (Canonical Datasets)",
           y = "Correlation with True HRF", x = "Dataset") +
      ylim(0, 1)
  }
  
  return(plots)
}


# Main execution
if (interactive()) {
  cat("Starting M-HRF-LSS benchmark evaluation...\n")
  cat("This will test the algorithm on multiple realistic datasets.\n\n")
  
  # Run benchmarks
  results <- run_all_benchmarks(presets = c("conservative", "balanced", "robust"))
  
  # Create visualizations
  plots <- plot_benchmark_summary(results)
  
  # Display plots
  print(plots$r2)
  print(plots$runtime)
  if (!is.null(plots$hrf)) print(plots$hrf)
  
  # Save results
  saveRDS(results, "benchmark_results.rds")
  cat("\nResults saved to benchmark_results.rds\n")
}