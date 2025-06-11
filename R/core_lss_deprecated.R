# Deprecated LSS functions - use core_lss.R instead

#' @title (DEPRECATED) Run LSS for Voxel Core
#' @description This function is deprecated. Use appropriate LS-A or LS-S functions instead.
#' @export
run_lss_for_voxel_core <- function(...) {
  .Deprecated(
    "run_lsa_for_voxel_corrected_full or run_lss_for_voxel",
    package = "manifold.hrf",
    old = "run_lss_for_voxel_core is deprecated. Use run_lsa_for_voxel_corrected_full() for LS-A or run_lss_for_voxel() for LS-S"
  )
  
  # Call the actual implementation
  stop("Function deprecated. Use run_lsa_for_voxel_corrected_full() or run_lss_for_voxel() instead.")
}

#' @title (DEPRECATED) Run LSS Voxel Loop Core
#' @description This function is deprecated. Use appropriate LS-A or LS-S functions instead.  
#' @export
run_lss_voxel_loop_core_deprecated <- function(...) {
  .Deprecated(
    "run_lsa_voxel_loop or run_lss_voxel_loop",
    package = "manifold.hrf",
    old = "run_lss_voxel_loop_core_deprecated is deprecated. Use run_lsa_voxel_loop() for LS-A or run_lss_voxel_loop() for LS-S"
  )
  
  # Call the actual implementation
  stop("Function deprecated. Use run_lsa_voxel_loop() or run_lss_voxel_loop() instead.")
}
