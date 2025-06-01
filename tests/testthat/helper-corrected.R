run_lss_voxel_loop_corrected_test <- function(Y_matrix,
                                             X_trial_onset_list,
                                             H_shapes_matrix,
                                             confounds_matrix,
                                             lambda = 1e-6) {
  P_conf <- prepare_projection_matrix(confounds_matrix, lambda)
  T_trials <- length(X_trial_onset_list)
  V <- ncol(Y_matrix)
  Beta <- matrix(0, T_trials, V)
  for (v in seq_len(V)) {
    y_proj <- as.vector(P_conf %*% Y_matrix[, v])
    Beta[, v] <- run_lss_for_voxel_corrected_full(
      Y_proj_voxel_vector = y_proj,
      X_trial_onset_list_of_matrices = X_trial_onset_list,
      H_shape_voxel_vector = H_shapes_matrix[, v],
      P_confound = P_conf,
      lambda_ridge = lambda
    )
  }
  Beta
}
