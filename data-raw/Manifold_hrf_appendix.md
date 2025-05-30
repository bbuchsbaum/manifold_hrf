Okay, this is a very well-structured and comprehensive proposal. I will now add an appendix detailing how its components would be implemented using the `neuroim2` and `fmrireg` R packages, as per your request.

---
**Revised Proposal: Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS) v3.0**

*(Original proposal content from Manifold.txt v3.0 as provided previously)...*
*(Sections 1 through 6 would be here)*
---

**Appendix: Implementation with `neuroim2` and `fmrireg` R Packages**

This appendix outlines how the M-HRF-LSS pipeline components can be implemented leveraging the functionalities of the `neuroim2` package for neuroimaging data structures and I/O, and the `fmrireg` package for HRF modeling, design matrix construction, and GLM fitting.

**A.1 General Package Roles:**

*   **`neuroim2`:**
    *   Data representation: `NeuroVol` for 3D volumes (masks, parameter maps), `NeuroVec` for 4D BOLD data or collections of volumes. `NeuroSpace` for spatial metadata. `LogicalNeuroVol` for masks.
    *   I/O: `read_vol()`, `read_vec()`, `write_vol()`, `write_vec()` for NIfTI files.
    *   Coordinate transformations: `index_to_grid()`, `grid_to_index()`, `coord_to_grid()`, `grid_to_coord()`.
    *   Basic image operations: arithmetic, subsetting, `series()` to extract voxel time series, `slices()` to extract 2D slices.
    *   Spatial operations: `gaussian_blur()`, `resample()`.
*   **`fmrireg`:**
    *   HRF modeling: `hrf()` function with various basis sets (e.g., `HRF_SPMG1`, `HRF_BSPLINE`), custom HRF creation with `as_hrf()`, HRF manipulation (`lag_hrf()`, `block_hrf()`), `evaluate.HRF()` for sampling HRFs.
    *   Design matrix construction: `event_model()` for event-related designs, `baseline_model()` for drift and nuisance regressors, `fmri_model()` to combine them, `design_matrix()` to get the matrix. `sampling_frame()` for run timing.
    *   GLM fitting: `fmri_lm()` for OLS estimation (voxel-wise or chunk-wise), `fmri_rlm()` for robust GLM. Provides coefficients, residuals, and statistics.
    *   LSS: Internal functions like `lss_fast()` support Least Squares Separate estimation, accessible via `estimate_betas(method="lss")`.
    *   Contrast estimation: `contrast_weights()`, `contrast_set()`.
    *   Simulation: `simulate_bold_signal()`.

**A.2 Component-wise Implementation Details:**

**Component -1: Pre-flight Data Quality Control (QC)**

*   **BOLD Data Input:**
    *   `neuroim2::read_vec(nifti_file_path)` loads the 4D BOLD data into a `NeuroVec` object.
    *   `neuroim2::read_vol(mask_file_path)` loads a brain mask into a `LogicalNeuroVol` or `NeuroVol`.
*   **TR & Timings:**
    *   TR is supplied as a scalar. `fmrireg::sampling_frame()` is constructed using TR and `run_length` (number of scans per run).
    *   Event timings are read from external files (e.g., CSV, text) into a `data.frame` to be used with `fmrireg::event_model()`.
*   **Motion Spikes & DVARS:**
    *   Motion parameters are typically external (e.g., from FSL MCFLIRT). Framewise displacement can be calculated from these.
    *   DVARS calculation would be custom R code operating on the `NeuroVec` data (e.g., using `neuroim2::series()` to get voxel data per run).
*   **Output:** `qc_flags` stored in a list.

**Component 0: HRF Manifold Construction**

*   **Inputs:**
    *   `L_library`: A standard R `matrix` (`p x N`), where each column is an HRF shape.
*   **Steps:**
    1.  **Affinity Matrix `W`:**
        *   Pairwise Euclidean distances: `stats::dist()` or `proxy::dist()` on `t(L_library)`.
        *   `k_local_nn` for `σ_i`: Loop through HRFs. For each HRF `i`, find its `k_local_nn`-th nearest neighbor using, e.g., `RANN::nn2(query=L_library[,i,drop=F], data=L_library, k=k_local_nn+1)$nn.dist[,k_local_nn+1]`.
        *   `W_ij`: Standard R matrix operations. If `N` is large, `W` can be a `Matrix::sparseMatrix`.
    2.  **Markov Matrix `S`:** Standard R matrix operations (`rowSums`, `diag`, `%*%`).
    3.  **Diffusion Map Coordinates `Φ_raw`:** `RSpectra::eigs_sym(S, k = m_manifold_dim_target + 5, which = "LM")`.
    4.  **Manifold Coordinates `Φ` & Automatic `m` selection:** Eigenvalue processing and selection in base R.
    5.  **HRF Reconstructor `B_reconstructor`:** Matrix algebra (`%*%`, `fmrireg::regularized_solve()` or base `solve(crossprod(Φ) + small_ridge * diag(ncol(Φ)))`).
*   **Outputs:** `B_reconstructor` (matrix), `Φ` (matrix), parameters.

**Component 1: Voxel-wise HRF Manifold Coordinate & Initial Condition Amplitude Estimation**

*   **Inputs:**
    *   `Y_bold`: `neuroim2::NeuroVec` object.
    *   `Z_confounds`: Optional matrix of nuisance regressors.
    *   `B_reconstructor`: Matrix from Component 0 (`p x m`).
    *   `h_ref_shape_canonical`: Vector, e.g., `fmrireg::HRF_SPMG1(seq(0, span, by=TR))`.
*   **Steps:**
    1.  **Confound Projection (Manual if needed, or include in GLM):**
        *   If manual: `Y_data <- as.matrix(series(Y_bold, mask_indices)); QZ <- qr.Q(qr(Z_confounds)); Y_proj_data <- Y_data - QZ %*% (t(QZ) %*% Y_data)`. Design matrices would also need projection.
        *   Alternatively, include `Z_confounds` in `fmrireg::baseline_model`.
    2.  **Per-Condition Design in Manifold Basis `Z_list` & Combined Design `X_tilde`:**
        *   The core idea `Xc %*% B_reconstructor` means convolving onsets of condition `c` with the `m` basis HRFs in `B_reconstructor`.
        *   Create a custom `fmrireg::HRF` object from `B_reconstructor`:
            *   `manifold_hrf_basis <- fmrireg::as_hrf(f = function(t_points) { apply(B_reconstructor, 2, function(basis_vec) { stats::spline(x = seq(0, (p-1)*TR_precision, by=TR_precision), y = basis_vec, xout = t_points, method="fmm")$y }) }, nbasis = m_manifold_dim, span = (p-1)*TR_precision)`
            *   Here `TR_precision` is the sampling rate of `L_library` HRFs.
        *   Construct `fmrireg::event_model` for all conditions using this `manifold_hrf_basis`:
            *   `event_data_df` needs a factor column `condition_factor` indicating which of the `k` conditions each event belongs to.
            *   `event_spec <- fmrireg::event_model(onset_col ~ hrf(condition_factor, basis = manifold_hrf_basis), data = event_data_df, block = ~run_col, sampling_frame = sframe)`
            *   `X_tilde_event_part <- fmrireg::design_matrix(event_spec)` (This will be `n x (km)`).
        *   If `Z_confounds` are handled via GLM:
            *   `baseline_spec <- fmrireg::baseline_model(nuisance_list = list_of_confound_matrices_per_run, sframe = sframe)`
            *   `fmri_model_spec <- fmrireg::fmri_model(event_spec, baseline_spec)`
            *   `X_combined_design <- fmrireg::design_matrix(fmri_model_spec)`
    3.  **GLM Solve for `Gamma_coeffs`:**
        *   `Y_masked_matrix <- neuroim2::series(Y_bold, which(mask_vol > 0))` (or `Y_proj_data` if projected).
        *   The `fmrireg::fmri_lm()` function can perform the voxel-wise GLM:
            *   `dset_for_lm <- fmrireg::matrix_dataset(Y_masked_matrix, TR=TR, run_length=run_lengths_vector, event_table=event_data_df_with_condition_factor)`
            *   `fit_results <- fmrireg::fmri_lm(model_formula = onset_col ~ hrf(condition_factor, basis = manifold_hrf_basis), block = ~run_col, baseline_model = baseline_spec_if_confounds_in_glm, dataset = dset_for_lm, strategy = "chunkwise", use_fast_path=TRUE, lambda_ridge = lambda_gamma)`
            *   `Gamma_coeffs_matrix <- fmrireg::coef(fit_results)`
            *   This matrix will be `(km + n_confounds) x V`. We need the `km x V` part corresponding to the event regressors. `fmrireg::term_indices(fit_results$model)` can identify these columns.
    4.  **Extract `ξ_raw_allvox` and `β_raw_allvox` via SVD (per voxel):**
        *   Loop `vx` from 1 to `V`:
            *   `G_vx_vector = Gamma_coeffs_event_part[, vx]`.
            *   `G_vx_matrix = matrix(G_vx_vector, nrow = m_manifold_dim, ncol = k)`.
            *   `svd_G_vx = base::svd(G_vx_matrix)`.
            *   Populate `Xi_raw_allvox` and `Beta_raw_allvox` (Base R matrix ops).
    5.  **Intrinsic Identifiability:**
        *   `xi_ref_coord = MASS::ginv(B_reconstructor) %*% h_ref_shape_canonical` (Base R).
        *   Loop through voxels, apply sign/scale adjustments using base R math.
*   **Outputs:** `Xi_ident_allvox` (matrix), `Beta_ident_allvox` (matrix). Store as `neuroim2::NeuroVec` if desired (each dimension/condition as a volume).

**Component 2: Spatial Smoothing of Manifold Coordinates**

*   **Inputs:**
    *   `Xi_ident_allvox`: `m x V` matrix.
    *   `voxel_coordinates_matrix`: `cds <- neuroim2::coords(mask_vol, real=FALSE)` gives `V x 3` grid coordinates.
*   **Steps:**
    1.  **`voxel_graph_laplacian_Lsp`:**
        *   k-NN search: `RANN::nn2(cds, cds, k=num_neighbors_Lsp+1)`.
        *   Build sparse adjacency matrix, then Laplacian: `Matrix::sparseMatrix()`, then `Diagonal(V, rowSums(adj)) - adj`.
    2.  **Smoothing Loop:** For each manifold dimension `j`:
        *   `b_vec = Xi_ident_allvox[j,]`.
        *   `A_mat = Matrix::Diagonal(V) + lambda_spatial_smooth * voxel_graph_laplacian_Lsp`.
        *   `Xi_smoothed_allvox[j,] = Matrix::solve(A_mat, b_vec)`.
*   **Outputs:** `Xi_smoothed_allvox` (matrix). Store as `neuroim2::NeuroVec`.

**Component 3: Trial-wise Amplitude Estimation (LSS)**

*   **Inputs:**
    *   `Y_proj`: `neuroim2::NeuroVec` (or original `Y_bold` if `A_lss_fixed` includes confounds).
    *   `H_shapes_allvox = B_reconstructor %*% Xi_smoothed_allvox` (`p x V` matrix).
    *   `A_lss_fixed`: Matrix of fixed LSS regressors (e.g., confounds + intercept).
*   **Steps (Voxel-Batched Woodbury - Implemented with custom loop):**
    1.  **Precompute LSS Fixed Components:** `P_lss`, `p_lss_vec` using base R math or `MASS::ginv`. Add jitter to `(AᵀA + λI)` if near-singular.
    2.  **Voxel Loop (parallelizable):** For each voxel `v`:
        a.  `hrf_v_shape <- H_shapes_allvox[,v]`.
        b.  Create `C_v` (`n x T` matrix): For each trial `t`, `C_v[,t]` is the regressor for trial `t` convolved with `hrf_v_shape`.
            *   This involves: `X_trial_onset_list[[t]] %*% hrf_v_shape`.
            *   `X_trial_onset_list[[t]]` is the `n x p` matrix of `p` delta functions for trial `t`.
            *   Alternatively, and more `fmrireg`-like: define `hrf_obj_v <- fmrireg::as_hrf(function(time_pts) { stats::spline(x=seq(0,(p-1)*TR_p,by=TR_p), y=hrf_v_shape, xout=time_pts)$y }, nbasis=1, span=(p-1)*TR_p)`.
            *   `trial_event_model_v <- fmrireg::event_model(onset_col ~ trialwise(basis=hrf_obj_v), data=event_table_with_trial_indices, block=~run_col, sampling_frame=sframe)`.
            *   `C_v <- fmrireg::design_matrix(trial_event_model_v)`. This `C_v` would be `n x T`.
        c.  **Woodbury LSS Core:** Apply the matrix algebra from proposal section "Component 3, Step 5b" to `Y_proj_voxel_v = neuroim2::series(Y_proj, v_index)`, `C_v`, `A_lss_fixed`, `P_lss`, `p_lss_vec`. This yields `Beta_trial_allvox_for_voxel_v`.
*   **Outputs:** `Beta_trial_allvox` (`T x V` matrix). Store as `neuroim2::NeuroVec`.
    *   *Note:* `fmrireg::estimate_betas(method="lss")` could be used if `A_lss_fixed` is empty or just an intercept, and if its internal LSS matches "Mode B". If `A_lss_fixed` is complex, a custom loop implementing the Woodbury steps is more direct.

**Component 4: Alternating Optimization & Final Condition Beta Re-estimation**

*   **Inputs:** `Y_proj` (`NeuroVec`), `H_shapes_allvox` (`p x V` matrix).
*   **Steps (Voxel Loop, default 1 pass):** For each voxel `v`:
    a.  `hrf_v_shape <- H_shapes_allvox[,v]`.
    b.  `hrf_obj_v <- fmrireg::as_hrf(...)` as in Component 3.
    c.  `event_model_final_v <- fmrireg::event_model(onset_col ~ hrf(condition_factor, basis=hrf_obj_v), data=event_data_df, block=~run_col, sampling_frame=sframe)`.
    d.  `fmri_model_final_v <- fmrireg::fmri_model(event_model_final_v, baseline_model_for_confounds_if_any)`.
    e.  `Y_voxel_data_matrix <- matrix(neuroim2::series(Y_proj, v_index), ncol=1)`.
    f.   `dset_voxel <- fmrireg::matrix_dataset(Y_voxel_data_matrix, TR=TR, run_length=run_lengths, event_table=event_data_df)`.
    g.  `fit_voxel <- fmrireg::fmri_lm(model=fmri_model_final_v, dataset=dset_voxel, lambda_ridge = lambda_beta_final)`.
    h.  `Beta_condition_final_allvox[,v] <- fmrireg::coef(fit_voxel)[condition_indices]`.
*   **Outputs:** `Beta_condition_final_allvox` (`k x V` matrix). Store as `neuroim2::NeuroVec`.

**A.3 Validation Framework (VAL-SIM-01 Integration)**

*   **Synthetic BOLD Data Generation:**
    *   Define true HRF shapes using `fmrireg`'s predefined HRFs (e.g., `HRF_SPMG1`, `HRF_GAUSSIAN`) or custom `as_hrf()` objects.
    *   Use `fmrireg::simulate_bold_signal()`: this function takes onsets, condition labels, amplitudes per trial, a list of HRFs (one per condition, or one for all), TR, and total time to generate `Y_clean`.
        *   `fmrireg::simulate_bold_signal(onsets_vector, condition_labels_vector, amplitudes_matrix_TxC, list_of_true_hrfs, TR, total_time)`
    *   Add noise: `Y_noisy_matrix = Y_clean_matrix + matrix(rnorm(n*V, mean=0, sd=noise_sd_value), n, V)`.
*   **Running Pipeline:** Use the `matrix_dataset` class in `fmrireg` to wrap the synthetic `Y_noisy_matrix` and feed it into the M-HRF-LSS pipeline.
*   **Metrics:** RMSE, correlation, Dice coefficient calculated using base R functions on the true vs. estimated beta maps/HRF parameters.

**A.4 Pipeline Quality Control (QC) Report (QC-REPORT-01 Integration)**

*   **Data Access for Plotting:**
    *   Manifold coordinates (`Xi_smoothed_allvox`), HRF shapes (`H_reconstructed_smoothed_allvox`), Betas (`Beta_condition_final_allvox`, `Beta_trial_allvox`) can be converted to `neuroim2::NeuroVol` or `NeuroVec` objects for easy slicing and mapping.
    *   `neuroim2::slices(neuro_vol_object, along_axis_index, slice_indices)` to get 2D slices.
    *   `neuroim2::plot(neuro_slice_object)` or `neuroim2::plot(neuro_vol_object)` (which plots multiple slices).
*   **Statistical Info:** R² values can be extracted from `fmri_lm` fit objects or calculated manually.
*   **Report Generation:** Use `rmarkdown::render()` with a template that incorporates plots and summary statistics.

**A.5 Data Structures:**

*   **BOLD Data:** `neuroim2::NeuroVec` (4D). `neuroim2::NeuroSpace` stores geometry.
*   **Masks:** `neuroim2::LogicalNeuroVol` or `neuroim2::NeuroVol` (3D).
*   **Parameter Maps (ξ, β):** `neuroim2::NeuroVec` where each "volume" is a map for one manifold dimension or one condition/trial.
*   **HRF Shapes (`L_library`, `B_reconstructor`, `H_shapes_allvox`):** Standard R `matrix`.
*   **Design Matrices (`X_tilde`, `C_v`, `A_lss_fixed`):** Standard R `matrix` or `Matrix::sparseMatrix`. `fmrireg` produces design matrices as `tbl_df` or `matrix`.
*   **fMRI Model Specifications:** `fmrireg` uses S3 classes like `fmri_model`, `event_model`, `baseline_model`, `hrfspec`.

**A.6 Computational Considerations:**

*   **Voxel-wise loops:** `fmrireg`'s `fmri_lm` with `strategy="chunkwise"` handles voxel batching and can utilize `future.apply` for parallelization if a `future` plan is set. Custom loops for SVD (Component 1) or LSS (Component 3) can be parallelized using `future.apply::future_lapply` or `parallel::mclapply`.
*   **Sparse Matrices:** `Matrix` package for graph Laplacians.
*   **BLAS/LAPACK:** R's matrix operations are typically linked against optimized BLAS/LAPACK libraries.

This appendix provides a roadmap for implementing the M-HRF-LSS pipeline by mapping its algorithmic steps to the functionalities provided by the `neuroim2` and `fmrireg` packages, supplemented by base R capabilities and specialized packages like `RSpectra` and `Matrix`.