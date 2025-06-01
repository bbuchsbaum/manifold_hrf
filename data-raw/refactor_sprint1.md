Okay, I've reviewed the `manifoldhrf` proposal (v3.0) and the provided `fmrireg` codebase. The core idea is to leverage `fmrireg`'s robust HRF handling and event modeling capabilities to feed into the `manifoldhrf` pipeline, making `manifoldhrf` more compatible and potentially more powerful.

Here's a breakdown of integration points and proposed tickets, structured for an engineer with limited `fmrireg` knowledge, including small code snippets where helpful.

**Overall Integration Strategy:**

1.  **HRF Library Definition & Sampling:** `manifoldhrf` will primarily use `fmrireg`'s HRF objects and `fmrireg::evaluate()` to generate the `L_library_matrix` for manifold construction.
2.  **"Raw" Event Matrix Generation:** `manifoldhrf` needs `n x p` (time x HRF samples) design matrices for each condition/trial *before* HRF convolution by the manifold method. We need to explore how `fmrireg::event_model` can produce these (e.g., by specifying an identity/FIR-like basis to `fmrireg`).
3.  **Output Compatibility:** Estimated HRFs from `manifoldhrf` should ideally be representable as `fmrireg` HRF objects for easier use in subsequent `fmrireg` analyses or plotting.

---

**Proposed Tickets for Integration:**

---

**Ticket MHRF-FMRIREG-001: Enhance `construct_hrf_manifold_nim` to Accept `fmrireg` HRF Objects**

*   **Objective:** Modify `manifoldhrf/R/neuroimaging_wrappers.R::construct_hrf_manifold_nim` to accept a list of `fmrireg` HRF objects (e.g., `HRF_SPMG1`, `HRF_GAMMA`, custom HRFs created with `fmrireg::as_hrf`) as `hrf_library_source`.
*   **Current State (manifoldhrf):** `construct_hrf_manifold_nim` has stubs for "FLOBS", "half_cosine", "gamma_grid" and can take a raw matrix. It doesn't directly process `fmrireg` HRF objects.
    ```R
    # In manifoldhrf/R/neuroimaging_wrappers.R (simplified)
    # L_library_matrix <- switch(hrf_library_source,
    #   "FLOBS" = create_flobs_library(...),
    #   "gamma_grid" = create_gamma_grid_library(...),
    #   ...
    # )
    # if (is.list(hrf_library_source)) { ... L_library_matrix <- matrix(0, p, N); ... }
    ```
*   **Desired State (fmrireg integration):**
    *   If `hrf_library_source` is a list, check if elements are `fmrireg::HRF` objects.
    *   If they are, use `fmrireg::evaluate()` to sample each HRF at `TR_precision` over `hrf_duration` to build `L_library_matrix`.
*   **Guidance (fmrireg):**
    *   `fmrireg` HRF objects are functions with S3 class "HRF".
    *   Use `inherits(obj, "HRF")` to check.
    *   Use `fmrireg::evaluate(hrf_obj, time_points_vector)` to get the shape.
    ```R
    # Snippet for engineer:
    # if (is.list(hrf_library_source) && all(sapply(hrf_library_source, inherits, "HRF"))) {
    #   time_points <- seq(0, hrf_duration, by = TR_precision)
    #   L_library_list <- lapply(hrf_library_source, function(hrf_fmrireg) {
    #     fmrireg::evaluate(hrf_fmrireg, time_points)
    #   })
    #   L_library_matrix <- do.call(cbind, L_library_list)
    #   # Ensure normalization if required by manifoldhrf (e.g., L_library_matrix <- apply(L_library_matrix, 2, function(h) h / sum(abs(h))) )
    # }
    ```
*   **Files to Update:**
    *   `manifoldhrf/R/neuroimaging_wrappers.R`

---

**Ticket MHRF-FMRIREG-002: Adapt `manifoldhrf` internal HRF library generators to output `fmrireg` HRF objects**

*   **Objective:** Modify internal library generators in `manifoldhrf` (e.g., `create_gamma_grid_library`, `create_flobs_library`) to return a list of `fmrireg::HRF` objects instead of directly returning a matrix. This promotes consistency. `construct_hrf_manifold_nim` would then always use `fmrireg::evaluate` (as per MHRF-FMRIREG-001) on these generated objects.
*   **Current State (manifoldhrf):** Functions like `create_gamma_grid_library` in `manifoldhrf/R/neuroimaging_wrappers.R` directly compute and return a `p x N` matrix.
*   **Desired State (fmrireg integration):** These functions should generate individual HRF functions and wrap them using `fmrireg::as_hrf()`.
*   **Guidance (fmrireg):**
    *   Use `fmrireg::as_hrf(fun, name, nbasis, span, params)` to create HRF objects.
    ```R
    # Snippet for engineer (e.g., modifying create_gamma_grid_library):
    # Original: L[, i] <- dgamma(time_points, shape = grid$shape[i], scale = grid$scale[i])
    # New:
    # hrf_func_i <- function(t) { stats::dgamma(t, shape = grid$shape[i], scale = grid$scale[i]) }
    # hrf_obj_i <- fmrireg::as_hrf(hrf_func_i,
    #                              name = paste0("gamma_shape", grid$shape[i], "_scale", grid$scale[i]),
    #                              span = hrf_duration)
    # library_list[[i]] <- hrf_obj_i
    # ...
    # return(library_list) # To be consumed by construct_hrf_manifold_nim
    ```
*   **Files to Update:**
    *   `manifoldhrf/R/neuroimaging_wrappers.R`

---

**Ticket MHRF-FMRIREG-003: Represent `B_reconstructor` as an `fmrireg` Multi-Basis HRF Object**

*   **Objective:** The output `B_reconstructor` from `manifoldhrf` Component 0 (Manifold Construction) defines a new basis set. Package this `p x m` matrix as an `fmrireg` compatible multi-basis HRF object.
*   **Current State (manifoldhrf):** `B_reconstructor` is a raw matrix. `manifoldhrf/R/neuroimaging_wrappers.R::create_manifold_hrf_object` is a placeholder.
*   **Desired State (fmrireg integration):** `create_manifold_hrf_object` should use `fmrireg::as_hrf` or `fmrireg::bind_basis` to represent `B_reconstructor`. Each column of `B_reconstructor` can be treated as a basis function.
*   **Guidance (fmrireg):**
    *   `fmrireg::as_hrf()` can take a function. We need to create a function that, given `t` and a basis index `j`, returns the `j`-th column of `B_reconstructor` evaluated/interpolated at `t`. More simply, if `B_reconstructor` columns are already sampled HRF shapes, `fmrireg::bind_basis` can be used on `fmrireg::empirical_hrf` objects.
    ```R
    # Snippet for engineer (inside create_manifold_hrf_object):
    # Assume B_reconstructor is p_samples x m_dims
    # time_points_for_B <- seq(0, by = TR_precision, length.out = nrow(B_reconstructor))
    #
    # basis_functions_list <- list()
    # for (j in 1:ncol(B_reconstructor)) {
    #   # Create an empirical HRF for each column of B_reconstructor
    #   basis_j_shape <- B_reconstructor[, j]
    #   # fmrireg::empirical_hrf takes time and values
    #   basis_functions_list[[j]] <- fmrireg::empirical_hrf(
    #                                   time_points_for_B,
    #                                   basis_j_shape,
    #                                   name = paste0(name, "_basis", j)
    #                                 )
    # }
    # manifold_hrf_basis <- do.call(fmrireg::bind_basis, basis_functions_list)
    # attr(manifold_hrf_basis, "name") <- name # Set overall name
    # return(manifold_hrf_basis)
    ```
    *   The `nbasis` attribute should be `m_manifold_dim`.
*   **Files to Update:**
    *   `manifoldhrf/R/neuroimaging_wrappers.R`

---

**Ticket MHRF-FMRIREG-004: Generate `n x p` "Raw" Event Design Matrices using `fmrireg`**

*   **Objective:** `manifoldhrf` (Components 1, 3, 4) requires `n x p` design matrices for conditions (`X_condition_list`) and trials (`X_trial_list`), where `p` is the HRF length (number of samples). These matrices essentially represent time-shifted delta functions (or boxcars for non-zero durations) *before* convolution with an estimated HRF. This ticket is to replace manual construction of these matrices with an `fmrireg`-idiomatic approach.
*   **Current State (manifoldhrf):**
    *   `manifoldhrf/R/mhrf_lss.R::.create_design_matrices` manually creates these by shifting `diag(actual_length)`.
    *   Tests like `test-alternating-optimization.R` and `test-fmrireg-benchmarks.R` also show manual construction.
*   **Desired State (fmrireg integration):**
    *   Utilize `fmrireg::event_model` with a special "identity" or "FIR-like" basis that produces the desired `n x p` matrix structure for each event/condition.
    *   The `p` (HRF length) parameter from `manifoldhrf` needs to map to the number of basis functions in this FIR-like `fmrireg` HRF.
*   **Guidance (fmrireg):**
    *   `fmrireg` has `HRF_FIR` defined in `R/hrf-functions.R`. This is the key.
        ```R
        # HRF_FIR generator:
        # hrf_fir_generator <- function(nbasis = 12, span = 24)
        # It creates 'nbasis' boxcar functions, each lasting 'span/nbasis'.
        # We need 'p' basis functions, each effectively a delta function for one TR slice of the HRF.
        # So, nbasis should be `p` (manifoldhrf's HRF length in samples).
        # The span should be `p * TR_of_HRF_sampling` (e.g., `p * TR_precision` if HRFs are finely sampled).
        # If TR_of_HRF_sampling is same as BOLD TR, then span = `p * TR`.
        ```
    *   Modify `.create_design_matrices` in `manifoldhrf/R/mhrf_lss.R`:
        1.  Define an `fmrireg::HRF` object that represents the `p` delta functions. This could be `HRF_FIR(nbasis = p, span = p * TR_precision_for_HRF_output)`. *Crucially, `TR_precision_for_HRF_output` must be the TR at which `manifoldhrf` expects its `p`-sample HRFs to be defined.* If this is the same as the acquisition TR, then `span = p * TR`.
        2.  Use this FIR-like HRF within an `fmrireg::event_model` formula.
        3.  Extract the resulting design matrices. `fmrireg::term_matrices(event_model_object)` might give the required `n x p` matrices if each "condition" in `manifoldhrf` corresponds to a term. If `X_condition_list` in `manifoldhrf` is for *actual experimental conditions*, and *each* needs an `n x p` matrix, then the `fmrireg` formula would be like `onset ~ hrf(experimental_condition, basis=my_fir_basis)`. `fmrireg` would then produce an `n x (num_levels * p)` matrix. This needs to be split correctly.
        4.  Alternatively, `fmrireg::convolve.event_term` is used if `hrfspec` is attached. The goal is for `manifoldhrf` to pass an *identity* `hrfspec` here for this stage.

    ```R
    # Snippet for engineer (conceptual, inside .create_design_matrices):
    # p_hrf_samples <- params$p_hrf # HRF length in samples (e.g., 25)
    # TR_for_HRF_output <- TR # Assuming manifoldhrf p_hrf is at acquisition TR
                             # OR: TR_for_HRF_output <- params$TR_precision (if p_hrf is at finer scale)

    # Create the FIR-like basis that generates p columns, each representing a time point of the HRF
    # Each "basis" function of HRF_FIR is a boxcar. We want delta-like.
    # A simpler way for manifoldhrf's needs might be to use fmrireg::regressor with a simple identity HRF *manually* for each event onset.
    # However, to get an N x P matrix *per condition*, we'd use fmrireg:

    # Option A: If fmrireg can produce the n x p matrix *per condition level* directly.
    #           This is the ideal scenario.
    # my_fir_basis <- fmrireg::HRF_FIR(nbasis = p_hrf_samples, span = p_hrf_samples * TR_for_HRF_output)
    #
    # Example for one experimental condition 'cond_A_events' (data frame with onset, duration, condition='A')
    # sframe_obj <- fmrireg::sampling_frame(blocklens=n_timepoints, TR=TR) # for one run
    # ev_model_for_cond_A <- fmrireg::event_model(onset ~ hrf(condition, basis=my_fir_basis),
    #                                          data=cond_A_events, sampling_frame=sframe_obj)
    # X_cond_A_matrix <- fmrireg::design_matrix(ev_model_for_cond_A)
    # This X_cond_A_matrix would be n_timepoints x p_hrf_samples. This is what manifoldhrf needs.
    # This needs to be done for each *original* condition that manifoldhrf expects in X_condition_list.

    # Option B: If fmrireg gives N x (num_levels * P), then split it.
    # my_fir_basis <- fmrireg::HRF_FIR(nbasis = p_hrf_samples, span = p_hrf_samples * TR_for_HRF_output)
    # full_event_model <- fmrireg::event_model(onset ~ hrf(original_condition_factor, basis=my_fir_basis), data=events_df, ...)
    # X_full_fir <- fmrireg::design_matrix(full_event_model)
    # Then, X_condition_list would be derived by splitting X_full_fir columns.
    # e.g., if 'original_condition_factor' has levels L1, L2:
    # X_condition_list[[1]] <- X_full_fir[, 1:p_hrf_samples] (for L1)
    # X_condition_list[[2]] <- X_full_fir[, (p_hrf_samples+1):(2*p_hrf_samples)] (for L2)

    # This needs careful handling of how `manifoldhrf` defines "conditions" vs. `fmrireg` terms.
    # The `create_trial_matrices` function in `mhrf_lss_interface.R` seems to be a good place to adapt.
    # It currently uses `fmrireg::evaluate(trial_reg, block_times - block_times[1])`.
    # If `trial_reg`'s HRF is set to the `p`-sample FIR basis, this might work.
    ```
*   **Challenge:** `fmrireg`'s `design_matrix()` typically returns *convolved* regressors. We need `manifoldhrf` to do its *own* HRF application (using `B_reconstructor`). So, the "basis" passed to `fmrireg::hrf()` must be one that effectively creates the `n x p` Toeplitz-like matrix. `HRF_FIR(nbasis=p, span=p*TR)` is the most promising candidate if each "basis" of `HRF_FIR` corresponds to one time-shifted delta (or narrow boxcar).
*   **Files to Update:**
    *   `manifoldhrf/R/mhrf_lss.R` (specifically `.create_design_matrices`)
    *   Potentially `manifoldhrf/R/mhrf_lss_interface.R` (`create_trial_matrices`)
    *   Relevant test files in `manifoldhrf/tests/testthat/`
See Appendix for a practical approach.

---

**Ticket MHRF-FMRIREG-005: Review and Refine `mhrf_lss` User Interface for HRF Library Specification**

*   **Objective:** Ensure the main `mhrf_lss` function (in `manifoldhrf/R/mhrf_lss_interface.R`, and the newer `mhrf_analyze` in `manifoldhrf/R/mhrf_lss.R`) correctly handles `hrf_library` arguments that can now be `fmrireg` HRF objects or lists thereof, passing them appropriately to the updated `construct_hrf_manifold_nim`.
*   **Current State (manifoldhrf):**
    *   `mhrf_lss_interface.R::mhrf_lss` takes `hrf_library` as string ("canonical", "flobs", "gamma_grid") or matrix.
    *   `mhrf_lss.R::mhrf_analyze` takes `hrf_library` as string ("auto", "spmg1", "gamma", "custom matrix").
*   **Desired State (fmrireg integration):**
    *   Both interfaces (or the preferred one, `mhrf_analyze`) should consistently accept:
        *   Predefined strings (e.g., "spmg1", "gamma_grid" which map to `fmrireg` or internal generators now returning `fmrireg` objects).
        *   A single `fmrireg::HRF` object.
        *   A list of `fmrireg::HRF` objects.
        *   A raw `p x N` matrix (legacy support).
    *   The logic within `mhrf_lss` or `.prepare_hrf_library` should correctly dispatch to `construct_hrf_manifold_nim` (or its core components).
*   **Guidance (fmrireg):** This is mostly about `manifoldhrf` internal consistency after MHRF-FMRIREG-001 and MHRF-FMRIREG-002 are implemented.
*   **Files to Update:**
    *   `manifoldhrf/R/mhrf_lss_interface.R` (if still primary)
    *   `manifoldhrf/R/mhrf_lss.R` (likely `.prepare_hrf_library` helper)

---

**Ticket MHRF-FMRIREG-006: Package Estimated Voxel-wise HRFs as `fmrireg` Empirical HRFs**

*   **Objective:** The final HRF estimates from `manifoldhrf` (`H_reconstructed_smoothed_allvox`, a `p x V` matrix) should be convertible into a list of `V` `fmrireg::empirical_hrf` objects, or a custom multi-HRF object for `fmrireg`.
*   **Current State (manifoldhrf):** Output is a raw matrix.
*   **Desired State (fmrireg integration):** Provide a helper function or method to convert the `p x V` matrix of HRF shapes into `fmrireg` compatible HRF objects.
*   **Guidance (fmrireg):**
    *   `fmrireg::empirical_hrf(time_points, hrf_values, name)` creates an HRF object from sampled data.
    ```R
    # Snippet for engineer:
    # Assume H_shapes_allvox is p x V, TR_precision_for_HRF_output is known
    # time_points_for_hrf <- seq(0, by = TR_precision_for_HRF_output, length.out = nrow(H_shapes_allvox))
    #
    # estimated_fmrireg_hrfs <- list()
    # for (v_idx in 1:ncol(H_shapes_allvox)) {
    #   estimated_fmrireg_hrfs[[v_idx]] <- fmrireg::empirical_hrf(
    #                                           time_points_for_hrf,
    #                                           H_shapes_allvox[, v_idx],
    #                                           name = paste0("voxel_", v_idx, "_mhrf")
    #                                         )
    # }
    # This list could be an output field or used to create a custom fmrireg basis.
    ```
*   **Files to Update:**
    *   Likely in `manifoldhrf/R/mhrf_result_methods.R` or a new utility file. Could be a method for the `mhrf_result` object.

---

**General Guidance for the Engineer:**

*   **Familiarize with `fmrireg` HRF System:** Understand `HRF` S3 class, `as_hrf`, `gen_hrf`, `bind_basis`, `evaluate.HRF`, `empirical_hrf`. Key files: `fmrireg/R/hrf.R`, `fmrireg/R/hrf-functions.R`, `fmrireg/R/hrf_decorators.R`.
*   **Familiarize with `fmrireg` Event Modeling:** Understand `event_model`, `design_matrix.event_model`, `term_matrices`. Key files: `fmrireg/R/event_model.R`, `fmrireg/R/event_vector.R` (for how terms are built).
*   **Focus on `manifoldhrf/R/neuroimaging_wrappers.R::construct_hrf_manifold_nim`** and **`manifoldhrf/R/mhrf_lss.R::.create_design_matrices`** as primary integration points for HRF library input and raw event matrix generation, respectively.
*   The `manifoldhrf/R/mhrf_lss_interface.R` file seems to be an older interface; the newer `mhrf_analyze` in `manifoldhrf/R/mhrf_lss.R` might be the target for user-facing changes. Clarify which is the canonical interface to update.
*   Pay close attention to sampling rates: `manifoldhrf`'s `TR_precision` for HRF sampling vs. `fmrireg`'s `TR` for BOLD acquisition. Ensure consistency when evaluating HRFs and constructing design matrices.
*   The `manifoldhrf` tests (e.g., `test-fmrireg-benchmarks.R`, `test-neuroimaging-wrappers.R`, `test-mhrf-lss-interface.R`) will need significant updates to reflect these `fmrireg` integrations.

This set of tickets should provide a structured approach to integrating `manifoldhrf` with `fmrireg`. The most challenging part will likely be MHRF-FMRIREG-004 due to the different ways `manifoldhrf` and `fmrireg` typically handle the step from event onsets to design matrices used for HRF application/estimation. But see Appendix for a practical approach.


Appendix:

Here's how we can structure this:

**1. Define `HRF_RAW_EVENT_BASIS` within `manifoldhrf`:**

This function will be placed in `manifoldhrf`, likely in a utility file or in `R/neuroimaging_wrappers.R` since it's related to interfacing with `fmrireg`-style HRF objects.

```R
# In manifoldhrf (e.g., R/neuroimaging_wrappers.R or a new R/fmrireg_helpers.R)

#' Create an HRF Basis for Raw Time Series Extraction (Temporary Location)
#'
#' @description
#' IMPORTANT: This function is temporarily located in the `manifoldhrf` package.
#' It is intended to be moved to the `fmrireg` package in a future update.
#'
#' This function creates an `fmrireg::HRF` object that, when used with
#' `fmrireg::hrf()` in an `fmrireg::event_model` formula, generates a design
#' matrix suitable for representing a raw, `p_length`-sample time course
#' aligned to event onsets. Each of the `p_length` basis functions
#' corresponds to one sample of this raw time course.
#'
#' @param p_length Integer, the number of samples in the desired raw HRF.
#' @param TR_sample Numeric, the time resolution of these `p_length` samples.
#'   This defines the duration of each individual "boxcar" basis function.
#' @param name Optional character string for the name of the HRF object.
#'   If `NULL`, a default name is generated.
#'
#' @return An `fmrireg::HRF` object.
#' @keywords internal
#' # Add @export if you want manifoldhrf users to directly call manifoldhrf::HRF_RAW_EVENT_BASIS
#' # For now, keeping internal as it's a helper for a specific pipeline step.
HRF_RAW_EVENT_BASIS <- function(p_length, TR_sample, name = NULL) {
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    stop("The 'fmrireg' package is required to create this HRF object.", call. = FALSE)
  }
  if (!is.numeric(p_length) || length(p_length) != 1 || p_length < 1 || p_length != round(p_length)) {
    stop("'p_length' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(TR_sample) || length(TR_sample) != 1 || TR_sample <= 0) {
    stop("'TR_sample' must be a positive numeric value.", call. = FALSE)
  }

  nbasis_for_fir <- as.integer(p_length)
  span_for_fir <- nbasis_for_fir * TR_sample
  bin_width <- TR_sample # Each "basis" of FIR is one TR_sample wide

  # Define the internal FIR-like function
  # This function will be wrapped by fmrireg::as_hrf
  # t_eval_relative_to_onset: time points relative to event onset
  f_fir_internal <- function(t_eval_relative_to_onset) {
    if (!is.numeric(t_eval_relative_to_onset) || length(t_eval_relative_to_onset) == 0) {
      return(matrix(0, nrow = 0, ncol = nbasis_for_fir))
    }
    output_matrix <- matrix(0, nrow = length(t_eval_relative_to_onset), ncol = nbasis_for_fir)

    for (i in seq_along(t_eval_relative_to_onset)) {
      current_t <- t_eval_relative_to_onset[i]
      # Check if current_t is within the span [0, span_for_fir)
      if (!is.na(current_t) && current_t >= 0 && current_t < span_for_fir) {
        # Determine which bin current_t falls into. bin_index is 1-based.
        # If current_t is exactly 0, it falls into the first bin.
        # If bin_width is 0, this will lead to issues, but TR_sample > 0 is checked.
        bin_index <- if (current_t == 0 && bin_width > 0) 1 else floor(current_t / bin_width) + 1
        # Cap bin_index at nbasis_for_fir to handle t_eval_relative_to_onset that might be exactly span_for_fir
        bin_index <- min(bin_index, nbasis_for_fir)
        output_matrix[i, bin_index] <- 1
      }
      # Values outside the span [0, span_for_fir) remain 0
    }
    return(output_matrix)
  }

  base_name <- paste0("RawEventBasis_p", p_length, "_tr", gsub("\\.", "p", as.character(TR_sample)))
  final_name <- name %||% base_name # Use manifoldhrf's %||%

  # Create the fmrireg::HRF object
  raw_event_hrf_obj <- fmrireg::as_hrf(
    f = f_fir_internal,
    name = final_name,
    nbasis = nbasis_for_fir,
    span = span_for_fir,
    params = list(p_length = p_length, TR_sample = TR_sample,
                  original_fir_nbasis = nbasis_for_fir,
                  original_fir_span = span_for_fir,
                  original_fir_bin_width = bin_width)
  )
  return(raw_event_hrf_obj)
}

# Ensure manifoldhrf's %||% is defined (e.g., in R/utils.R)
# `%||%` <- function(a, b) if (is.null(a)) b else a
```

**2. How `manifoldhrf` would use this (revised for clarity):**

The usage within `manifoldhrf`'s `.create_design_matrices` (or similar logic in `mhrf_lss_interface.R`'s `create_trial_matrices`) would be:

```R
# Inside manifoldhrf's .create_design_matrices or equivalent:

# Parameters available within manifoldhrf:
# p_hrf_samples: manifoldhrf's desired HRF length in samples (its 'p', e.g., 25 from params$p_hrf)
# TR_for_manifold_hrf: The TR at which these 'p' samples are defined.
#                      This could be the acquisition TR (params$TR) or a finer
#                      TR_precision (e.g., params$TR_precision for oversampled HRFs).
#                      Let's assume it's params$TR for simplicity here.
# acquisition_TR: The actual TR of the BOLD data (params$TR).
# events_df_with_condition_column: The event table for manifoldhrf.
# condition_factor_name: The column name in events_df_with_condition_column that defines conditions.
# run_column_name: The column name for run/block identifiers.
# n_timepoints_per_run_vector: Vector of scan counts per run (e.g., dataset$run_length).

# 1. Create the raw event basis HRF object using the function now defined in manifoldhrf
#    This HRF_RAW_EVENT_BASIS is now manifoldhrf::HRF_RAW_EVENT_BASIS
my_raw_hrf_object <- HRF_RAW_EVENT_BASIS(p_length = p_hrf_samples,
                                         TR_sample = TR_for_manifold_hrf) # TR_sample is key here

# 2. Prepare for fmrireg::event_model.
#    The `basis` argument to `fmrireg::hrf()` can take an actual HRF object.
sframe_obj <- fmrireg::sampling_frame(blocklens = n_timepoints_per_run_vector, TR = acquisition_TR)

# Construct the hrfspec to pass to fmrireg::event_model via its list interface
# (or construct a formula and ensure my_raw_hrf_object is in the eval env)
# Using the list interface is cleaner for passing objects.

# Get the symbolic representation of the condition factor
condition_sym <- rlang::sym(condition_factor_name)

# Create an hrfspec object using fmrireg::hrf()
# The `basis` argument takes our custom HRF object.
hrf_spec_for_raw_designs <- fmrireg::hrf(!!condition_sym, basis = my_raw_hrf_object)
# This `hrf_spec_for_raw_designs` is now an `hrfspec` object.

# Create the event model using the list interface for clarity
model_spec_raw <- fmrireg::event_model(
  formula_or_list = list(raw_events = hrf_spec_for_raw_designs), # Give the term a name
  data = events_df_with_condition_column,
  block = rlang::new_formula(NULL, rlang::sym(run_column_name)), # e.g., ~ run
  sampling_frame = sframe_obj
)

# 3. Extract the full design matrix produced by fmrireg
#    This matrix will be n_total_scans x (num_condition_levels * p_hrf_samples)
X_full_raw <- fmrireg::design_matrix(model_spec_raw)

# 4. Split X_full_raw into X_condition_list for manifoldhrf
#    The columns in X_full_raw are named by fmrireg based on the term tag ("raw_events"),
#    condition tags (e.g., "condition.LevelA"), and basis suffixes ("_b01", "_b02", ...).
#    We need to reliably map these back to manifoldhrf's conditions.

X_condition_list <- list()
# `fmrireg::conditions(model_spec_raw$terms$raw_events, expand_basis=FALSE)` gives base conditions
# `fmrireg::conditions(model_spec_raw$terms$raw_events, expand_basis=TRUE)` gives expanded names
# The columns of X_full_raw will correspond to the expanded names.

# Let's get the unique levels of the original condition factor
original_condition_levels <- levels(events_df_with_condition_column[[condition_factor_name]])
if (is.null(original_condition_levels)) original_condition_levels <- unique(as.character(events_df_with_condition_column[[condition_factor_name]]))


# The `conditions()` method for `event_term` (which is inside the `event_model`)
# combined with the `term_tag` and `basis_suffix` logic from `fmrireg`'s naming
# utility (`make_column_names`) is what generates the colnames of `X_full_raw`.
# We need to find which columns belong to each original condition level.

# Assuming fmrireg column names are like: "raw_events_condition.<LEVEL>_b<XX>"
# For each original_condition_level, select its p_hrf_samples columns.
term_tag_from_model <- names(model_spec_raw$terms)[1] # e.g., "raw_events"

for (level in original_condition_levels) {
  # Construct the expected prefix for this level's columns
  # fmrireg's `level_token` creates "FactorName.LevelName"
  condition_token <- fmrireg:::level_token(condition_factor_name, level) # Internal but shows intent
  # The full prefix before basis suffix:
  prefix_for_level <- paste0(term_tag_from_model, "_", condition_token)

  # Find columns matching this prefix
  # We need to ensure the grep is specific enough.
  # Columns are prefix_for_level_b01, prefix_for_level_b02, ... prefix_for_level_b<p_hrf_samples>
  level_cols_indices <- grep(paste0("^", prefix_for_level, "_b"), colnames(X_full_raw))

  if (length(level_cols_indices) == p_hrf_samples) {
    X_condition_list[[as.character(level)]] <- X_full_raw[, level_cols_indices, drop = FALSE]
  } else {
    warning(paste("Could not correctly extract columns for condition level:", level,
                  ". Expected", p_hrf_samples, "columns, found", length(level_cols_indices),
                  "matching prefix", prefix_for_level))
  }
}
# Make sure the list elements are named by the original condition levels
names(X_condition_list) <- original_condition_levels


# Now X_condition_list contains the n x p matrices manifoldhrf needs for each of its original conditions.
```

**Key changes and considerations:**

1.  **Location:** The `HRF_RAW_EVENT_BASIS` function is defined within `manifoldhrf`.
2.  **Dependency:** It explicitly calls `fmrireg::as_hrf`, so `fmrireg` must be an import.
3.  **FIR Logic:** The core logic of generating the FIR-like basis functions (`f_fir_internal`) is now part of this new function. This logic creates `p_length` distinct basis functions, each representing a "1" at a specific time lag relative to event onset, and "0" elsewhere within the HRF's sample window.
4.  **Naming:** `HRF_RAW_EVENT_BASIS` generates a descriptive name for the `fmrireg::HRF` object it creates.
5.  **Usage in `manifoldhrf`:**
    *   `manifoldhrf` first calls its *own* `HRF_RAW_EVENT_BASIS` to get an `fmrireg::HRF` object.
    *   This HRF object is then passed to `fmrireg::hrf(..., basis = my_raw_hrf_object)` when `manifoldhrf` sets up the `fmrireg::event_model`.
    *   `fmrireg` then uses this special HRF object to generate the `n x (num_levels * p_length)` design matrix.
    *   `manifoldhrf` will then need to correctly parse the column names of this matrix (which `fmrireg` generates based on its naming conventions: `term_tag_condition_tag_b##`) to split it into the `n x p_length` matrices for each of its original conditions. This splitting logic is crucial.

This approach isolates the specific FIR-like basis generation while still leveraging `fmrireg` for the main event modeling and design matrix construction. When `fmrireg` is updated, this temporary function in `manifoldhrf` can be removed, and `manifoldhrf` can call `fmrireg::HRF_RAW_EVENT_BASIS` directly.