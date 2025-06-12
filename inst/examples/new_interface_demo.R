# Demonstration of the new manifoldhrf interface
# This shows how to use the simplified API

library(manifoldhrf)

# Simulate some fMRI data
set.seed(123)
n_timepoints <- 200
n_voxels <- 100
n_trials <- 20
TR <- 2.0

# Create simulated BOLD data
fmri_data <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

# Create event information
events <- data.frame(
  onset = sort(runif(n_trials, 0, (n_timepoints - 20) * TR)),
  condition = rep(c("A", "B"), each = n_trials/2),
  duration = 0
)

# Example 1: Basic usage with default settings
# --------------------------------------------
cat("Example 1: Basic usage\n")
fit1 <- estimate_hrf_manifold(
  fmri_data = fmri_data,
  events = events,
  TR = TR,
  hrf_library = "fir"
)

# Print basic information
print(fit1)

# Get summary
summary(fit1)

# Extract coefficients
amplitudes <- coef(fit1, type = "amplitudes")
cat("\nCondition amplitudes (first 5 voxels):\n")
print(amplitudes[, 1:5])

# Example 2: Using different presets
# ----------------------------------
cat("\n\nExample 2: Fast estimation\n")

# Use fast preset for quick analysis
ctrl_fast <- manifold_control(preset = "fast")
print(ctrl_fast)

fit2 <- estimate_hrf_manifold(
  fmri_data = fmri_data,
  events = events,
  TR = TR,
  hrf_library = "fir",
  control = ctrl_fast
)

# Example 3: Custom control parameters
# ------------------------------------
cat("\n\nExample 3: Custom parameters\n")

# Create custom control with specific settings
ctrl_custom <- manifold_control(
  preset = "balanced",
  m_manifold_dim = 5,             # Use 5 manifold dimensions
  lambda_spatial_smooth = 1.0,     # Strong spatial smoothing
  max_iter = 200,                  # More iterations
  verbose_level = 2,               # Detailed output
  generate_qc_report = TRUE        # Generate QC metrics
)

fit3 <- estimate_hrf_manifold(
  fmri_data = fmri_data,
  events = events,
  TR = TR,
  hrf_library = "fir",
  control = ctrl_custom
)

# Example 4: Working with results
# -------------------------------
cat("\n\nExample 4: Working with results\n")

# Extract different components
trial_amps <- coef(fit3, type = "trial_amplitudes")
manifold_coords <- coef(fit3, type = "manifold_coords")

cat(sprintf("Trial amplitudes shape: %d x %d\n", nrow(trial_amps), ncol(trial_amps)))
cat(sprintf("Manifold coordinates shape: %d x %d\n", nrow(manifold_coords), ncol(manifold_coords)))

# Get fitted values and residuals
fitted_values <- fitted(fit3)
residuals <- residuals(fit3)

# Calculate R-squared
ss_tot <- sum((fmri_data - mean(fmri_data))^2)
ss_res <- sum(residuals^2)
r_squared <- 1 - ss_res/ss_tot
cat(sprintf("Overall R-squared: %.3f\n", r_squared))

# Example 5: Plotting results (if available)
# -----------------------------------------
cat("\n\nExample 5: Visualization\n")

# Plot HRF shapes for first few voxels
if (interactive()) {
  plot(fit3, type = "hrf", voxels = 1:3)
  
  # Plot manifold coordinates
  plot(fit3, type = "manifold")
  
  # Plot amplitude distribution
  plot(fit3, type = "amplitudes")
}

# Example 6: Preset comparison
# ----------------------------
cat("\n\nExample 6: Comparing presets\n")

presets <- c("fast", "balanced", "thorough")
for (preset in presets) {
  ctrl <- manifold_control(preset = preset)
  cat(sprintf("\n%s preset:\n", preset))
  cat(sprintf("  Max iterations: %d\n", ctrl$max_iter))
  cat(sprintf("  Manifold neighbors: %d\n", ctrl$num_neighbors_mfd))
  cat(sprintf("  Spatial neighbors: %d\n", ctrl$num_neighbors_Lsp))
  cat(sprintf("  Lambda (manifold): %.2f\n", ctrl$lambda_manifold))
  cat(sprintf("  Lambda (spatial): %.2f\n", ctrl$lambda_spatial_smooth))
}

cat("\n\nDemo complete!\n")