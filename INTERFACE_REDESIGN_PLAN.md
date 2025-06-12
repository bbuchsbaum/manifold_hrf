# Interface Redesign Plan for manifoldhrf

Based on discussion with Gemini, this document outlines the plan to create a unified, user-friendly interface for the manifoldhrf package.

## Overview

Transform the complex manifoldhrf package into a modern R package with a clean, consistent interface that follows best practices while retaining all advanced functionality.

## Core Design Decisions

### 1. Primary Function: `estimate_hrf_manifold()`

**Signature:**
```r
estimate_hrf_manifold(
  fmri_data,    # NeuroVec/array/matrix
  events,       # data.frame of events
  TR,           # numeric repetition time
  hrf_library,  # string preset or custom matrix
  mask = NULL,  # optional brain mask
  control = manifold_control()  # algorithm parameters
)
```

### 2. Control Function Pattern

Implement `manifold_control()` as a proper function with validation:

```r
manifold_control <- function(
  preset = "balanced",
  # Manifold parameters
  m_manifold_dim = NULL,
  lambda_manifold = 0.1,
  num_neighbors_mfd = 75,
  
  # Spatial smoothing
  lambda_spatial_smooth = 0.5,
  num_neighbors_Lsp = 6,
  
  # Optimization
  max_iter = 100,
  tolerance = 1e-6,
  
  # Computational
  n_cores = 1,
  verbose_level = 1,
  ...
) {
  # Apply presets
  params <- switch(preset,
    "fast" = list(max_iter = 50, num_neighbors_mfd = 50, tolerance = 1e-4),
    "balanced" = list(max_iter = 100, num_neighbors_mfd = 75, tolerance = 1e-6),
    "thorough" = list(max_iter = 500, num_neighbors_mfd = 100, tolerance = 1e-8),
    stop("Unknown preset: ", preset)
  )
  
  # Override with user-specified values
  user_params <- list(...)
  params[names(user_params)] <- user_params
  
  # Validation
  if (params$max_iter <= 0) stop("max_iter must be positive")
  if (params$lambda_spatial_smooth < 0) stop("lambda_spatial_smooth must be non-negative")
  
  structure(params, class = "manifold_control")
}
```

### 3. Output Object: `manifold_hrf_fit`

```r
structure(
  list(
    # Primary results
    amplitudes = beta_condition_final,        # Condition-level amplitudes
    trial_amplitudes = beta_trial,           # Trial-wise amplitudes
    hrf_shapes = H_shapes,                   # Estimated HRF shapes
    
    # Model fit
    fitted_values = fitted,                  # Predicted BOLD signal
    residuals = Y - fitted,                  # Residuals
    
    # Model details (nested for cleanliness)
    model_specific = list(
      manifold_coords = Xi_final,            # Manifold coordinates
      manifold = manifold_obj,               # Manifold object
      amplitudes_initial = beta_initial,     # Initial estimates
      spatial_laplacian = L_spatial,         # Spatial graph
      convergence_info = conv_info           # Iterations, tolerance achieved
    ),
    
    # Metadata
    call = match.call(),                     # For reproducibility
    control = control_params_used,           # Actual parameters used
    data_info = list(
      n_timepoints = n,
      n_voxels = V,
      n_trials = T,
      n_conditions = K,
      TR = TR
    ),
    
    # Optional components
    qc_metrics = qc_results                  # If QC was requested
  ),
  class = "manifold_hrf_fit"
)
```

## S3 Methods Implementation

### Essential Methods

```r
# 1. Print method - concise summary
print.manifold_hrf_fit <- function(x, ...) {
  cat("Manifold HRF Fit\n")
  cat("----------------\n")
  cat("Call: ", deparse(x$call), "\n")
  cat(sprintf("Data: %d timepoints, %d voxels, %d trials\n", 
              x$data_info$n_timepoints, x$data_info$n_voxels, x$data_info$n_trials))
  cat(sprintf("Conditions: %d\n", x$data_info$n_conditions))
  cat(sprintf("Manifold dimension: %d\n", ncol(x$model_specific$manifold_coords)))
  cat("\nCondition amplitudes:\n")
  print(summary(x$amplitudes))
  invisible(x)
}

# 2. Summary method - detailed information
summary.manifold_hrf_fit <- function(object, ...) {
  structure(
    list(
      call = object$call,
      data_info = object$data_info,
      amplitudes_summary = summary(object$amplitudes),
      trial_amplitude_summary = summary(object$trial_amplitudes),
      convergence = object$model_specific$convergence_info,
      r_squared = 1 - sum(object$residuals^2) / sum((object$fitted_values + object$residuals)^2)
    ),
    class = "summary.manifold_hrf_fit"
  )
}

# 3. Coefficient extraction
coef.manifold_hrf_fit <- function(object, type = c("amplitudes", "trial_amplitudes", "manifold_coords"), ...) {
  type <- match.arg(type)
  switch(type,
    amplitudes = object$amplitudes,
    trial_amplitudes = object$trial_amplitudes,
    manifold_coords = object$model_specific$manifold_coords
  )
}

# 4. Plotting methods
plot.manifold_hrf_fit <- function(x, type = c("hrf", "manifold", "amplitudes"), ...) {
  type <- match.arg(type)
  switch(type,
    hrf = plot_hrf_shapes(x$hrf_shapes, ...),
    manifold = plot_manifold_coords(x$model_specific$manifold_coords, ...),
    amplitudes = plot_amplitude_distribution(x$amplitudes, x$trial_amplitudes, ...)
  )
}

# 5. Prediction method
predict.manifold_hrf_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted_values)
  }
  # Generate predicted BOLD signal for new events
  # This requires convolving new events with estimated HRFs
  predict_bold_from_events(newdata, object$hrf_shapes, object$amplitudes, object$data_info$TR)
}

# 6. Standard accessors
fitted.manifold_hrf_fit <- function(object, ...) object$fitted_values
residuals.manifold_hrf_fit <- function(object, ...) object$residuals
```

### Advanced Methods (for tidyverse integration)

```r
# For broom compatibility
tidy.manifold_hrf_fit <- function(x, ...) {
  data.frame(
    condition = names(x$amplitudes),
    amplitude = x$amplitudes,
    n_trials = table(extract_conditions(x)),
    stringsAsFactors = FALSE
  )
}

glance.manifold_hrf_fit <- function(x, ...) {
  data.frame(
    n_timepoints = x$data_info$n_timepoints,
    n_voxels = x$data_info$n_voxels,
    n_trials = x$data_info$n_trials,
    n_conditions = x$data_info$n_conditions,
    manifold_dim = ncol(x$model_specific$manifold_coords),
    r_squared = calculate_r_squared(x),
    converged = x$model_specific$convergence_info$converged,
    iterations = x$model_specific$convergence_info$iterations
  )
}
```

## Transition Strategy

### Phase 1: Implement New Interface (v0.5.0)
1. Create `estimate_hrf_manifold()` function
2. Implement `manifold_control()` with validation
3. Create `manifold_hrf_fit` class and all S3 methods
4. Write comprehensive tests for new interface

### Phase 2: Deprecation Wrapper (v0.5.0)
```r
#' @export
mhrf_analyze <- function(...) {
  lifecycle::deprecate_warn(
    "0.5.0",
    "mhrf_analyze()",
    "estimate_hrf_manifold()",
    details = "See ?estimate_hrf_manifold for the new interface."
  )
  
  # Map old arguments to new
  args <- list(...)
  new_args <- translate_old_args(args)
  
  # Call new function
  result <- do.call(estimate_hrf_manifold, new_args)
  
  # Convert to old format for compatibility
  convert_to_mhrf_result(result)
}
```

### Phase 3: Documentation Update
1. Update all vignettes to use new interface
2. Create migration guide vignette
3. Update README with new examples
4. Add deprecation notices to old function docs

### Phase 4: Soft Deprecation (v0.6.0)
- Continue supporting old interface with warnings
- All new features only in new interface

### Phase 5: Hard Deprecation (v1.0.0)
- Remove old interface entirely
- Clean break for major version

## Implementation Checklist

- [ ] Create `R/estimate_hrf_manifold.R`
- [ ] Create `R/manifold_control.R`
- [ ] Create `R/manifold_hrf_fit_class.R` with class definition
- [ ] Implement all S3 methods in `R/manifold_hrf_fit_methods.R`
- [ ] Create transition wrapper in `R/mhrf_analyze_deprecated.R`
- [ ] Write tests in `tests/testthat/test-new-interface.R`
- [ ] Update NAMESPACE exports
- [ ] Create migration vignette
- [ ] Update all examples in documentation
- [ ] Add lifecycle badges to documentation
- [ ] Update NEWS.md with detailed change log

## Benefits

1. **User-Friendly**: Clear function names, intuitive parameter organization
2. **Extensible**: Easy to add new control parameters without breaking changes
3. **Standard**: Follows R modeling object conventions
4. **Modern**: Compatible with tidyverse workflows
5. **Maintainable**: Clear separation of concerns, easier to test

## Timeline

- Month 1: Implement core interface and methods
- Month 2: Create comprehensive tests and documentation
- Month 3: Beta release with deprecation warnings
- Month 4-6: Gather feedback and refine
- Month 7: Release v0.5.0 with new interface