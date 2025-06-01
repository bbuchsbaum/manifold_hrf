Okay, this is a great directive. We'll structure the sprint plan to clearly delineate between the "core computational engine" (which will operate on standard R matrices, lists, and `Matrix` sparse matrices) and the "neuroimaging layer" (which will handle `neuroim2`/`fmrireg` specific data structures, NIfTI I/O, and higher-level neuroimaging logic).

This promotes modularity, testability of the core algorithms with synthetic non-neuroimaging data, and flexibility for potential future backends or interfaces.

Here's a revised sprint plan reflecting this:

---

**Sprint Plan: M-HRF-LSS Pipeline Implementation (MVP - Core Engine First)**

**Sprint Goal:**
1.  Deliver a functional R implementation of the **core computational engine** for M-HRF-LSS (Components 0-4), operating on generic R data structures, focusing on algorithmic correctness and modularity.
2.  Provide initial wrapper functions demonstrating how to use this core engine with `neuroim2` and `fmrireg` data structures for a basic end-to-end run.

**User Stories (Guiding the Sprint):**
*   As a methods developer, I need robust, testable core functions for each M-HRF-LSS component that accept and return standard R data types (matrix, list, sparseMatrix).
*   As an fMRI methods developer, I need a thin neuroimaging wrapper layer that can prepare data from `neuroim2`/`fmrireg` objects, pass it to the core engine, and convert results back into neuroimaging formats.
*   As an fMRI researcher (early adopter), I want to run the wrapped M-HRF-LSS pipeline on a sample dataset.

---

**Tickets:**

**EPIC 0: Core - HRF Manifold Construction Engine (Component 0)**

1.  **`MHRF-CORE-MANIFOLD-01`: HRF Library Affinity & Markov Matrix (Core)**
    *   **Task:** `calculate_manifold_affinity_core(L_library_matrix, k_local_nn_for_sigma, use_sparse_W_params)`
        *   Input: `L_library_matrix` (`p x N` R matrix), `k_local_nn_for_sigma` (int), `use_sparse_W_params` (list for N > threshold, k_nn for sparse W).
        *   Output: `S_markov_matrix` (R matrix or `Matrix::sparseMatrix`).
        *   Details: Pairwise Euclidean distances, self-tuning bandwidth (`σ_i`, `σ_j`), `W_ij = exp(-dists_ij² / (σ_i * σ_j))`, optional k-NN sparsification for `W`, `S = D_inv %*% W`.
    *   **DoD:** Function returns `S_markov_matrix`. Tested with small matrix `L_library_matrix`.

2.  **`MHRF-CORE-MANIFOLD-02`: Diffusion Map Basis & Reconstructor (Core)**
    *   **Task:** `get_manifold_basis_reconstructor_core(S_markov_matrix, L_library_matrix, m_manifold_dim_target, m_manifold_dim_min_variance)`
        *   Input: `S_markov_matrix` (`N x N`), `L_library_matrix` (`p x N`), `m_manifold_dim_target` (int), `m_manifold_dim_min_variance` (numeric, e.g., 0.95).
        *   Output: List containing `B_reconstructor_matrix` (`p x m_final`), `Phi_coords_matrix` (`N x m_final`), `eigenvalues_S_vector`, `m_final_dim` (int), `m_auto_selected_dim` (int).
        *   Details: `RSpectra::eigs_sym` for eigenvectors, automatic `m` selection based on variance explained, `B_reconstructor = L_library %*% Φ %*% solve(crossprod(Φ) + ridge)`.
    *   **DoD:** Function returns list with specified matrices and scalars. Tested.

**EPIC 1: Core - Voxel-wise Manifold Fit Engine (Component 1)**

3.  **`MHRF-CORE-VOXFIT-01`: Confound Projection (Core)**
    *   **Task:** `project_out_confounds_core(Y_data_matrix, X_list_of_matrices, Z_confounds_matrix)`
        *   Input: `Y_data_matrix` (`n x V`), `X_list_of_matrices` (list of `k` `n x p` matrices), `Z_confounds_matrix` (`n x q_confound` or `NULL`).
        *   Output: List containing `Y_proj_matrix`, `X_list_proj_matrices`.
        *   Details: Uses `qr.Q(qr(Z_confounds_matrix))` for projection.
    *   **DoD:** Module tested for correct projection on matrices.

4.  **`MHRF-CORE-VOXFIT-02`: Design in Manifold Basis (Core)**
    *   **Task:** `transform_designs_to_manifold_basis_core(X_condition_list_proj_matrices, B_reconstructor_matrix)`
        *   Input: `X_condition_list_proj_matrices` (list of `k` `n x p` matrices), `B_reconstructor_matrix` (`p x m`).
        *   Output: `Z_list_of_matrices` (list of `k` `n x m` matrices).
    *   **DoD:** Module tested.

5.  **`MHRF-CORE-VOXFIT-03`: Main GLM for Gamma Coefficients (Core)**
    *   **Task:** `solve_glm_for_gamma_core(Z_list_of_matrices, Y_proj_matrix, lambda_gamma, orthogonal_approx_flag)`
        *   Input: `Z_list_of_matrices` (list of `k` `n x m` matrices), `Y_proj_matrix` (`n x V`), `lambda_gamma` (scalar), `orthogonal_approx_flag` (bool).
        *   Output: `Gamma_coeffs_matrix` (`(km) x V` R matrix).
        *   Details: Form `X_tilde` from `Z_list`, `XtX_tilde_reg = crossprod(X_tilde) + lambda_gamma * diag(k*m)`, solve for `Gamma_coeffs`. Implement optional `orthogonal_approx_flag` modification.
    *   **DoD:** `Gamma_coeffs_matrix` computed correctly.

6.  **`MHRF-CORE-VOXFIT-04`: SVD-based Extraction of Raw Manifold Coords & Amplitudes (Core)**
    *   **Task:** `extract_xi_beta_raw_svd_core(Gamma_coeffs_matrix, m_manifold_dim, k_conditions)`
        *   Input: `Gamma_coeffs_matrix` (`(km) x V`), `m_manifold_dim` (int), `k_conditions` (int).
        *   Output: List containing `Xi_raw_matrix` (`m x V`), `Beta_raw_matrix` (`k x V`).
        *   Details: Per-column SVD of reshaped `Gamma_coeffs`. Robustness for near-zero singular values.
    *   **DoD:** Raw `Xi_raw_matrix` and `Beta_raw_matrix` extracted.

7.  **`MHRF-CORE-VOXFIT-05`: Intrinsic Identifiability (Core)**
    *   **Task:** `apply_intrinsic_identifiability_core(Xi_raw_matrix, Beta_raw_matrix, B_reconstructor_matrix, h_ref_shape_vector, ident_scale_method, ident_sign_method)`
        *   Input: `Xi_raw_matrix` (`m x V`), `Beta_raw_matrix` (`k x V`), `B_reconstructor_matrix` (`p x m`), `h_ref_shape_vector` (`p x 1`), `ident_scale_method` (char), `ident_sign_method` (char).
        *   Output: List containing `Xi_ident_matrix` (`m x V`), `Beta_ident_matrix` (`k x V`).
        *   Details: `xi_ref_coord = ginv(B_reconstructor) %*% h_ref_shape_canonical`. Sign/scale adjustments per column.
    *   **DoD:** Identifiability constraints correctly applied to matrices.

**EPIC 2: Core - Spatial Smoothing Engine (Component 2)**

8.  **`MHRF-CORE-SPSMOOTH-01`: Graph Laplacian Construction (Core)**
    *   **Task:** `make_voxel_graph_laplacian_core(voxel_coords_matrix, num_neighbors_Lsp)`
        *   Input: `voxel_coords_matrix` (`V x 3` R matrix), `num_neighbors_Lsp` (int).
        *   Output: `L_sp_sparse_matrix` (`V x V` `Matrix::sparseMatrix`).
        *   Details: Use `RANN::nn2` for k-NN search on coordinates, construct sparse adjacency, then Laplacian.
    *   **DoD:** Sparse `L_sp_sparse_matrix` generated.

9.  **`MHRF-CORE-SPSMOOTH-02`: Apply Spatial Smoothing to Manifold Coordinates (Core)**
    *   **Task:** `apply_spatial_smoothing_core(Xi_ident_matrix, L_sp_sparse_matrix, lambda_spatial_smooth)`
        *   Input: `Xi_ident_matrix` (`m x V`), `L_sp_sparse_matrix` (`V x V`), `lambda_spatial_smooth` (scalar).
        *   Output: `Xi_smoothed_matrix` (`m x V` R matrix).
        *   Details: Loop `m` times, solve `(I_V + λ_sp L_sp) ξ_j_smooth = ξ_j_ident` using `Matrix::solve`.
    *   **DoD:** `Xi_smoothed_matrix` computed.

**EPIC 3: Core - Trial-wise LSS Engine (Component 3)**

10. **`MHRF-CORE-LSS-01`: LSS Fixed Regressor Precomputation (Core)**
    *   **Task:** `prepare_lss_fixed_components_core(A_lss_fixed_matrix, intercept_col_index_in_Alss, lambda_ridge_Alss)`
        *   Input: `A_lss_fixed_matrix` (`n x q_lss`), `intercept_col_index_in_Alss` (int or `NULL`), `lambda_ridge_Alss` (scalar).
        *   Output: List containing `P_lss_matrix` (`q_lss x n`), `p_lss_vector` (`n x 1`).
        *   Details: `P_lss = solve(AᵀA + λI) Aᵀ`. Handle `p_lss_vec` based on intercept. Add jitter for `AᵀA` solve.
    *   **DoD:** `P_lss_matrix` and `p_lss_vector` returned.

11. **`MHRF-CORE-LSS-02`: Voxel-Specific HRF Shape Reconstruction (Core)**
    *   **Task:** `reconstruct_hrf_shapes_core(B_reconstructor_matrix, Xi_smoothed_matrix)`
        *   Input: `B_reconstructor_matrix` (`p x m`), `Xi_smoothed_matrix` (`m x V`).
        *   Output: `H_shapes_allvox_matrix` (`p x V` R matrix).
    *   **DoD:** Voxel-specific HRF shapes matrix computed.

12. **`MHRF-CORE-LSS-03`: Woodbury LSS Kernel for Single Voxel (Core)**
    *   **Task:** `run_lss_for_voxel_core(Y_proj_voxel_vector, X_trial_onset_list_of_matrices, H_shape_voxel_vector, A_lss_fixed_matrix, P_lss_matrix, p_lss_vector)`
        *   Input: `Y_proj_voxel_vector` (`n x 1`), `X_trial_onset_list_of_matrices` (list of `T` `n x p` matrices), `H_shape_voxel_vector` (`p x 1`), `A_lss_fixed_matrix` (`n x q_lss`), `P_lss_matrix`, `p_lss_vector`.
        *   Output: `beta_trial_voxel_vector` (`T x 1` R vector).
        *   Details: Forms `C_v` for the voxel, then Woodbury steps for `S_effective_regressors_v`, then `crossprod(S_effective_regressors_v, Y_proj_voxel_vector)`.
    *   **DoD:** Function returns `T x 1` trial betas for one voxel. Tested.

13. **`MHRF-CORE-LSS-04`: Main LSS Loop (Core)**
    *   **Task:** `run_lss_voxel_loop_core(Y_proj_matrix, X_trial_onset_list_of_matrices, H_shapes_allvox_matrix, A_lss_fixed_matrix, P_lss_matrix, p_lss_vector, ram_heuristic_GB_for_Rt)`
        *   Input: `Y_proj_matrix` (`n x V`), `X_trial_onset_list_of_matrices` (list of `T` `n x p`), `H_shapes_allvox_matrix` (`p x V`), `A_lss_fixed_matrix`, `P_lss_matrix`, `p_lss_vector`, `ram_heuristic_GB_for_Rt` (numeric).
        *   Output: `Beta_trial_allvox_matrix` (`T x V` R matrix).
        *   Details: Optionally precompute `R_t_allvox` matrices (`list` of `n x V` matrices) if RAM permits. Loop over voxels, calling `run_lss_for_voxel_core`.
    *   **DoD:** `Beta_trial_allvox_matrix` computed.

**EPIC 4: Core - Alternating Optimization & Final Betas (Component 4)**

14. **`MHRF-CORE-ALTOPT-01`: Final Condition Beta Re-estimation (Core)**
    *   **Task:** `estimate_final_condition_betas_core(Y_proj_matrix, X_condition_list_proj_matrices, H_shapes_allvox_matrix, lambda_beta_final, control_alt_list)`
        *   Input: `Y_proj_matrix` (`n x V`), `X_condition_list_proj_matrices` (list of `k` `n x p`), `H_shapes_allvox_matrix` (`p x V`), `lambda_beta_final` (scalar), `control_alt_list` (for iterations).
        *   Output: `Beta_condition_final_matrix` (`k x V` R matrix).
        *   Details: Implements Component 4, Step 2. For MVP, `max_iter=1`. Voxel-wise GLM using reconstructed HRFs.
    *   **DoD:** `Beta_condition_final_matrix` computed.

**EPIC 5: Neuroimaging Layer - Wrappers & I/O (MVP)**

15. **`MHRF-NIM-IO-MANIFOLD-01`: Enhanced Manifold Construction Wrapper (Neuroimaging Layer)**
    *   **Task:** `construct_hrf_manifold_nim(hrf_library_source, TR_precision, ...)`
        *   Input: `hrf_library_source` (e.g., path to file with HRFs, list of `fmrireg::HRF` objects, predefined library name like "FLOBS" or "half_cosine").
        *   Additional inputs: `TR_precision` (sampling rate for HRF evaluation), manifold parameters.
        *   Output: List with `B_reconstructor_matrix`, `manifold_hrf_basis` (custom `fmrireg::HRF` object), and other manifold params.
        *   **Details:**
            *   Support multiple HRF sources: `fmrireg::HRF_SPMG1`, `fmrireg::HRF_GAUSSIAN`, custom `fmrireg::as_hrf()` objects
            *   Convert `fmrireg::HRF` objects to `L_library_matrix` using `fmrireg::evaluate.HRF()`
            *   Create custom `manifold_hrf_basis` using `fmrireg::as_hrf()` from `B_reconstructor_matrix`
            *   Calls `MHRF-CORE-MANIFOLD-01` & `-02` with error handling and validation
    *   **DoD:** Wrapper produces manifold components and `fmrireg`-compatible HRF basis objects from neuroimaging-relevant sources.

16. **`MHRF-NIM-WRAP-SUBJECT-01`: Enhanced Subject-Level Processing Wrapper (Neuroimaging Layer)**
    *   **Task:** `process_subject_mhrf_lss_nim(bold_input, mask_input, event_input, confound_input, manifold_objects, params_list)`
        *   Input: 
            *   `bold_input` (NIfTI path, `neuroim2::NeuroVec`, or matrix)
            *   `mask_input` (path, `neuroim2::LogicalNeuroVol`, or logical vector)
            *   `event_input` (path to CSV/TSV, `data.frame` with onset/condition/run columns)
            *   `confound_input` (optional, path or matrix for nuisance regressors)
            *   `manifold_objects` (from `MHRF-NIM-IO-MANIFOLD-01`)
            *   `params_list` (all pipeline parameters)
        *   Output: List of R matrices plus processing metadata.
        *   **Enhanced Details:**
            *   **Data Loading:** `neuroim2::read_vec()` for BOLD, `neuroim2::read_vol()` for mask, handle multiple input formats
            *   **Coordinate Extraction:** `neuroim2::coords(mask_vol, real=FALSE)` for voxel coordinates (Component 2)
            *   **BOLD Matrix Preparation:** `neuroim2::series(Y_bold, which(mask_vol > 0))` to get `Y_data_matrix`
            *   **Design Matrix Construction:** 
                *   Create `fmrireg::sampling_frame()` from TR and run lengths
                *   Build `fmrireg::event_model()` using `manifold_hrf_basis` for conditions and trials
                *   Generate both condition-wise and trial-wise design matrices using `fmrireg::design_matrix()`
                *   Handle confound integration via `fmrireg::baseline_model()` if provided
            *   **Pipeline Orchestration:** Call core engine functions (`MHRF-CORE-...`) in correct sequence with proper error handling
            *   **Memory Management:** Implement RAM heuristics for large datasets, optional HDF5 backing
    *   **DoD:** Comprehensive wrapper handles multiple input formats, orchestrates full pipeline, includes robust error handling.

17. **`MHRF-NIM-OUTPUT-01`: Enhanced Results Packaging & Visualization (Neuroimaging Layer)**
    *   **Task:** Create comprehensive output object with `neuroim2` integration and visualization methods.
        *   Functions: `package_mhrf_results_nim()`, plus S3 methods for the result object.
        *   **Input:** `core_results_list`, `reference_space` (`neuroim2::NeuroSpace`), `mask_vol`, `original_inputs`, `processing_metadata`
        *   **Enhanced Output Object (`mhrf_results` S3 class):**
            *   **Core Matrices:** All pipeline outputs as matrices (for programmatic access)
            *   **NeuroImaging Objects:** Key results as `neuroim2::NeuroVol`/`NeuroVec` objects:
                *   `Xi_smoothed_volumes`: `neuroim2::NeuroVec` with one volume per manifold dimension
                *   `H_reconstructed_volumes`: `neuroim2::NeuroVec` with HRF time courses
                *   `Beta_condition_volumes`: `neuroim2::NeuroVec` with condition-wise activation maps
                *   `Beta_trial_volumes`: `neuroim2::NeuroVec` with trial-wise activation maps
            *   **Metadata:** Parameters, processing time, QC flags, software versions
        *   **S3 Methods:**
            *   `print.mhrf_results()`: Summary of key results and QC status
            *   `summary.mhrf_results()`: Detailed statistics and diagnostics
            *   `plot.mhrf_results()`: Key visualizations using `neuroim2::plot()` for volumes, base R for HRF shapes/eigenvalues
            *   `write_results.mhrf_results()`: Export results to NIfTI files
            *   `[.mhrf_results()`: Subsetting interface for accessing specific components
    *   **DoD:** Comprehensive output object with `neuroim2` integration, full S3 method suite, and export capabilities.

**EPIC 6: Validation, QC, and Documentation (MVP)**

18. **`MHRF-VALIDATE-SIM-01`: Synthetic Data Generation & Validation Script (Core & Neuroimaging)**
    *   **Task:** Create a script `run_mhrf_lss_simulation_core()`
        *   Generates synthetic `Y_matrix` (`n x V`), true HRF parameters (`m x V`), true trial betas (`T x V`), true condition betas (`k x V`) as R matrices/vectors.
        *   Uses `fmrireg::simulate_bold_signal()` to generate `Y_clean_matrix`, then add noise.
        *   Calls the *core* M-HRF-LSS pipeline functions with these matrices.
        *   Computes and reports metrics (RMSE, correlation) for HRF params, condition betas, trial betas.
    *   **DoD:** Simulation script runs, produces validation metrics against known ground truth matrices.

19. **`MHRF-QC-REPORT-01`: Comprehensive QC Report Generation (Neuroimaging Layer)**
    *   **Task:** R Markdown template for comprehensive HTML QC report.
        *   Inputs: Output object from `MHRF-NIM-OUTPUT-01`, pipeline parameters, and processing metadata.
        *   **Report Contents:**
            *   **Input Parameters & QC Flags Summary:** All pipeline parameters and any QC flags with red/yellow/green badges.
            *   **Manifold Diagnostics (DIAG-MANIFOLD-01):**
                *   Diffusion map eigenvalue spectrum (scree plot)
                *   Cumulative variance explained by `m` dimensions, chosen `m` vs `m_auto`
                *   Reconstruction error (`||L - BΦ||F / ||L||F`) vs `m` plot
                *   (Optional) Interactive 3D scatter plot of `Φ` if `m=3`
            *   **HRF & Beta Diagnostics:**
                *   Maps of `Xi_smoothed_allvox` (first 3 components) using `neuroim2::plot`
                *   Before/after spatial smoothing comparison maps of one `ξ` coordinate
                *   Histogram of HRF peak times and FWHM from `H_reconstructed_smoothed_allvox`
                *   Median voxel-wise R² of GLM fit (Component 1) and LSS fit (Component 3)
                *   Convergence plot for alternating optimization (ΔHRF/ΔBeta per iteration) if Component 4 `max_iter > 1`
            *   **QC Flags Summary:** Include flags for `low_trial_count`, `poor_fits` (e.g., mean `R2_vox < 0.1` for >30% voxels), `unstable_hrf` (e.g., HRF peaks < 2s or > 10s), `high_motion_spikes`, `high_noise_dvars`
            *   **Failure-Mode Checklist:** Guidance on interpreting flags and troubleshooting
            *   **Performance Metrics:** Elapsed time and iteration counts for key components
    *   **DoD:** Template generates comprehensive HTML report with all specified diagnostic plots and summary statistics.

20. **`MHRF-TEST-UNIT-01`: Unit Tests for Core Functions**
    *   **Task:** Write `testthat` tests for key core engine functions (e.g., `MHRF-CORE-MANIFOLD-01`, `MHRF-CORE-VOXFIT-03`, `MHRF-CORE-LSS-03`).
    *   Focus on input/output dimensions, algorithmic correctness with small, predictable matrices.
    *   **DoD:** Core functions have unit tests.

21. **`MHRF-DOC-MVP-01`: MVP Documentation & Example**
    *   **Task:** Document main wrapper functions (`construct_hrf_manifold`, `process_subject_mhrf_lss`) and the output object.
    *   Create a simple example script demonstrating the neuroimaging layer.
    *   **DoD:** Basic user documentation and a runnable example script.

---

**Sprint Review Focus:**
*   Correctness of core engine functions with matrix inputs.
*   Successful execution of the neuroimaging wrapper layer on a test dataset.
*   Plausibility of estimated HRFs, condition betas, and trial betas from the wrapped pipeline.
*   Initial performance benchmarks (CPU-based R).
*   Validation metrics from the simulation script.

**Post-MVP Considerations (Future Sprints):**
*   Full implementation of Component -1 (Pre-flight QC) with `neuroim2` and external tool wrappers.
*   Parallelization of core voxel loops (e.g., `MHRF-CORE-VOXFIT-04`, `MHRF-CORE-LSS-04`) using `future.apply`.
*   `RcppArmadillo` optimization for core computational bottlenecks if identified.
*   Advanced QC report elements and visualizations using `neuroim2` plotting.
*   Robust parameter selection helpers (e.g., for `lambda_spatial_smooth`, `m_manifold_dim`).

This separation clearly prioritizes the general computational engine, making it the primary deliverable of the MVP sprint, with the neuroimaging specifics layered on top for practical application.