run_lss_voxel_loop_corrected_test <- function(Y_matrix,
                                             X_trial_onset_list,
                                             H_shapes_matrix,
                                             confounds_matrix,
                                             lambda = 0) {
  # Use the new LS-A Woodbury implementation (was previously doing LS-S incorrectly)
  T_trials <- length(X_trial_onset_list)
  V <- ncol(Y_matrix)
  Beta <- matrix(0, T_trials, V)
  for (v in seq_len(V)) {
    Beta[, v] <- run_lsa_woodbury(
      Y_proj_voxel_vector = Y_matrix[, v],  # No pre-projection needed
      X_trial_onset_list_of_matrices = X_trial_onset_list,
      H_shape_voxel_vector = H_shapes_matrix[, v],
      Z_confounds = confounds_matrix,
      lambda = lambda
    )
  }
  Beta
}
