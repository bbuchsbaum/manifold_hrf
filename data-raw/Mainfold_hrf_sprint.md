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

15. **`MHRF-NIM-IO-01`: Manifold Construction Wrapper (Neuroimaging)**
    *   **Task:** `construct_hrf_manifold(hrf_library_source, m_dim_target, m_min_var, ...)`
        *   Input: `hrf_library_source` (e.g., path to file, list of `fmrireg::HRF` objects).
        *   Output: List containing `B_reconstructor_matrix` and other manifold parameters.
        *   Details: Loads/generates `L_library_matrix` from source, calls `MHRF-CORE-MANIFOLD-01` & `-02`. Handles caching of results.
    *   **DoD:** Wrapper produces and saves/caches manifold components.

16. **`MHRF-NIM-SUBJECT-PROCESS-01`: Subject-Level Processing Wrapper (Neuroimaging)**
    *   **Task:** `process_subject_mhrf_lss(bold_file, mask_file, event_file, confound_file, manifold_objects, params_list)`
        *   Input: File paths or `neuroim2`/`fmrireg` objects, precomputed `manifold_objects`, list of all pipeline parameters.
        *   Output: A structured list/object containing all key outputs (`Xi_raw_matrix`, `Beta_condition_initial_matrix`, `Xi_smoothed_matrix`, `H_reconstructed_smoothed_matrix`, `Beta_trial_matrix`, `Beta_condition_final_matrix`).
        *   Details:
            *   Loads BOLD (`NeuroVec`), mask (`LogicalNeuroVol`).
            *   Prepares `event_data_df`, `run_lengths`, `TR` for `fmrireg::sampling_frame`.
            *   Uses `fmrireg::event_model` to create initial design structures (pre-manifold projection).
            *   Converts `NeuroVec` data to `n x V` matrix using `neuroim2::series(Y_bold, which(mask_vol > 0))`.
            *   Converts `fmrireg` design matrices to simple R matrices to feed into core functions.
            *   Calls core engine functions for Components 1, 2, 3, 4.
            *   Handles interim data conversions between matrix and `NeuroVec` representations if needed (e.g., for `Xi_smoothed_matrix` to be visualized).
    *   **DoD:** End-to-end processing for one subject, producing all matrix-based outputs.

17. **`MHRF-NIM-OUTPUT-02`: Define and Implement M-HRF-LSS Output Object (Neuroimaging)**
    *   **Task:** Convert matrix outputs from `MHRF-NIM-SUBJECT-PROCESS-01` into `neuroim2::NeuroVol` or `neuroim2::NeuroVec` objects. Create an S3/S4 class to store these `neuroim2` objects and parameters.
    *   Implement basic `print()`, `summary()`, and `plot()` (e.g., plot mean reconstructed HRF as a `NeuroVol` map of one `ξ` coordinate).
    *   **DoD:** Structured output object with `neuroim2` data and basic utility methods.

**EPIC 6: Validation, QC, and Documentation (MVP)**

18. **`MHRF-VALIDATE-SIM-01`: Synthetic Data Generation & Validation Script (Core & Neuroimaging)**
    *   **Task:** Create a script `run_mhrf_lss_simulation_core()`
        *   Generates synthetic `Y_matrix` (`n x V`), true HRF parameters (`m x V`), true trial betas (`T x V`), true condition betas (`k x V`) as R matrices/vectors.
        *   Uses `fmrireg::simulate_bold_signal()` to generate `Y_clean_matrix`, then add noise.
        *   Calls the *core* M-HRF-LSS pipeline functions with these matrices.
        *   Computes and reports metrics (RMSE, correlation) for HRF params, condition betas, trial betas.
    *   **DoD:** Simulation script runs, produces validation metrics against known ground truth matrices.

19. **`MHRF-QC-REPORT-01`: Basic QC Report Elements (Neuroimaging)**
    *   **Task:** Draft R Markdown template.
    *   Include: Manifold eigenvalue spectrum plot, reconstruction error vs. `m` plot.
    *   Plot a few example smoothed HRF shapes (`H_reconstructed_smoothed_allvox`).
    *   Map of one `Xi_smoothed` component (as `NeuroVol`).
    *   QC flags from Component -1.
    *   **DoD:** Template generates an HTML report with specified plots from pipeline outputs.

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

This separation clearly prioritizes the general computational engine, making it the primary deliverable of the MVP sprint, with the neuroimaging specifics layered on top for practical application.Okay, this is an excellent way to structure the sprint plan. It focuses on building the core computational engine first, ensuring it's robust and testable with standard R data structures. The neuroimaging-specific parts are then layered on top. This approach is great for modularity and maintainability.

Here's the ticketed sprint list:

---

**Sprint Plan: M-HRF-LSS Pipeline Implementation (MVP - Core Engine First)**

**Sprint Goal:**
1.  Deliver a functional R implementation of the **core computational engine** for M-HRF-LSS (Components 0-4), operating on generic R data structures (R `matrix`, `list` of matrices, `Matrix::sparseMatrix`), focusing on algorithmic correctness and modularity.
2.  Provide initial wrapper functions demonstrating how to use this core engine with neuroimaging data, potentially using `neuroim2` and `fmrireg` for data handling and design matrix generation before passing matrices to the core.

**User Stories (Guiding the Sprint):**
*   As a methods developer, I need robust, testable core functions for each M-HRF-LSS component that accept and return standard R data types.
*   As an fMRI methods developer, I need a thin neuroimaging wrapper layer that can prepare data (e.g., from NIfTI files or existing R objects), pass it as matrices/lists to the core engine, and convert matrix results back into meaningful neuroimaging formats or summaries.
*   As an fMRI researcher (early adopter), I want to run the wrapped M-HRF-LSS pipeline on a sample dataset.

---

**Tickets:**

**EPIC 0: Core Engine - HRF Manifold Construction (Component 0)**

1.  **`MHRF-CORE-MANIFOLD-01`: HRF Library Affinity & Markov Matrix Calculation (Core)**
    *   **Task:** Create function `calculate_manifold_affinity_core(L_library_matrix, k_local_nn_for_sigma, sparse_W_threshold = 5000, k_nn_for_W_sparse = 20)`
        *   Input: `L_library_matrix` (`p x N` R matrix of HRF shapes).
        *   Output: `S_markov_matrix` (`N x N` R matrix, or `Matrix::sparseMatrix` if sparsified).
        *   Details: Compute pairwise Euclidean distances. Implement self-tuning bandwidth `σ_i`, `σ_j`. Compute affinity `W_ij = exp(-dists_ij² / (σ_i * σ_j))`. Optional k-NN sparsification for `W`. Compute Markov Matrix `S = D_inv %*% W`.
    *   **DoD:** Function returns `S_markov_matrix`. Tested with small `L_library_matrix`.

2.  **`MHRF-CORE-MANIFOLD-02`: Diffusion Map Eigendecomposition & Reconstructor `B` (Core)**
    *   **Task:** Create function `get_manifold_basis_reconstructor_core(S_markov_matrix, L_library_matrix, m_manifold_dim_target, m_manifold_dim_min_variance = 0.95)`
        *   Input: `S_markov_matrix`, `L_library_matrix`, `m_manifold_dim_target` (int), `m_manifold_dim_min_variance` (numeric).
        *   Output: List containing `B_reconstructor_matrix` (`p x m_final`), `Phi_coords_matrix` (`N x m_final`), `eigenvalues_S_vector`, `m_final_dim` (int), `m_auto_selected_dim` (int). Ensure consistent eigenvector signs (e.g. first element positive).
        *   Details: Use `RSpectra::eigs_sym`. Automatic `m` selection. `B_reconstructor = L_library_matrix %*% Phi_coords_matrix %*% solve(crossprod(Phi_coords_matrix) + 1e-8 * diag(m_final_dim))`.
    *   **DoD:** Function returns list with matrices and scalars. Tested.

**EPIC 1: Core Engine - Voxel-wise Manifold Fit (Component 1)**

3.  **`MHRF-CORE-VOXFIT-01`: Confound Projection (Core)**
    *   **Task:** Create function `project_out_confounds_core(Y_data_matrix, X_list_of_matrices, Z_confounds_matrix = NULL)`
        *   Input: `Y_data_matrix` (`n x V` R matrix), `X_list_of_matrices` (list of `k` `n x p` R matrices), `Z_confounds_matrix` (optional `n x q_confound` R matrix).
        *   Output: List containing `Y_proj_matrix` and `X_list_proj_matrices`.
        *   Details: Uses `qr.Q(qr(Z_confounds_matrix, LAPACK=TRUE))` for projection operator `Q_Z`. `Y_proj = Y - Q_Z %*% (t(Q_Z) %*% Y)`.
    *   **DoD:** Module tested for correct matrix projection.

4.  **`MHRF-CORE-VOXFIT-02`: Per-Condition Design in Manifold Basis (Core)**
    *   **Task:** Create function `transform_designs_to_manifold_basis_core(X_condition_list_proj_matrices, B_reconstructor_matrix)`
        *   Input: `X_condition_list_proj_matrices`, `B_reconstructor_matrix` (`p x m`).
        *   Output: `Z_list_of_matrices` (list of `k` `n x m` R matrices).
    *   **DoD:** Module tested.

5.  **`MHRF-CORE-VOXFIT-03`: Main GLM Solve for `Gamma_coeffs` (Core)**
    *   **Task:** Create function `solve_glm_for_gamma_core(Z_list_of_matrices, Y_proj_matrix, lambda_gamma, orthogonal_approx_flag = FALSE)`
        *   Input: `Z_list_of_matrices`, `Y_proj_matrix` (`n x V`), `lambda_gamma` (scalar), `orthogonal_approx_flag` (bool).
        *   Output: `Gamma_coeffs_matrix` (`(km) x V` R matrix).
        *   Details: Form `X_tilde = do.call(cbind, Z_list_of_matrices)`. Form `XtX_tilde_reg` with ridge penalty. Handle `orthogonal_approx_flag`. Solve `Gamma_coeffs = solve(XtX_tilde_reg, crossprod(X_tilde, Y_proj_matrix))`.
    *   **DoD:** `Gamma_coeffs_matrix` computed correctly.

6.  **`MHRF-CORE-VOXFIT-04`: SVD-based Extraction of `ξ_raw` and `β_raw` (Core)**
    *   **Task:** Create function `extract_xi_beta_raw_svd_core(Gamma_coeffs_matrix, m_manifold_dim, k_conditions)`
        *   Input: `Gamma_coeffs_matrix`, `m_manifold_dim` (int), `k_conditions` (int).
        *   Output: List with `Xi_raw_matrix` (`m x V`), `Beta_raw_matrix` (`k x V`).
        *   Details: Per-column SVD loop. Handle near-zero singular values.
    *   **DoD:** Raw `Xi` and `Beta` matrices extracted.

7.  **`MHRF-CORE-VOXFIT-05`: Intrinsic Identifiability for `ξ` and `β` (Core)**
    *   **Task:** Create function `apply_intrinsic_identifiability_core(Xi_raw_matrix, Beta_raw_matrix, B_reconstructor_matrix, h_ref_shape_vector, ident_scale_method = "l2_norm", ident_sign_method = "canonical_correlation")`
        *   Input: `Xi_raw_matrix`, `Beta_raw_matrix`, `B_reconstructor_matrix`, `h_ref_shape_vector` (`p x 1` R vector for canonical HRF).
        *   Output: List with `Xi_ident_matrix` (`m x V`), `Beta_ident_matrix` (`k x V`).
        *   Details: Compute `xi_ref_coord` using `MASS::ginv(B_reconstructor_matrix) %*% h_ref_shape_vector`. Apply sign and scale adjustments per column.
    *   **DoD:** Identifiability constraints correctly applied to matrices.

**EPIC 2: Core Engine - Spatial Smoothing (Component 2)**

8.  **`MHRF-CORE-SPSMOOTH-01`: Graph Laplacian Construction (Core)**
    *   **Task:** Create function `make_voxel_graph_laplacian_core(voxel_coordinates_matrix, num_neighbors = 6)`
        *   Input: `voxel_coordinates_matrix` (`V x 3` R matrix of coordinates).
        *   Output: `voxel_graph_laplacian_Lsp` (`V x V` `Matrix::sparseMatrix`).
        *   Details: Use `RANN::nn2` for k-NN graph, then construct Laplacian.
    *   **DoD:** Sparse `L_sp` matrix generated.

9.  **`MHRF-CORE-SPSMOOTH-02`: Apply Spatial Smoothing to `ξ` Coordinates (Core)**
    *   **Task:** Create function `apply_spatial_smoothing_core(Xi_ident_matrix, voxel_graph_laplacian_Lsp, lambda_spatial_smooth)`
        *   Input: `Xi_ident_matrix` (`m x V`), `voxel_graph_laplacian_Lsp`, `lambda_spatial_smooth` (scalar).
        *   Output: `Xi_smoothed_matrix` (`m x V` R matrix).
        *   Details: Loop `m` times, solving sparse system `(Diagonal(V) + λ_sp L_sp) ξ_j_smooth = ξ_j_ident` using `Matrix::solve`.
    *   **DoD:** `Xi_smoothed_matrix` computed.

**EPIC 3: Core Engine - Trial-wise LSS (Component 3)**

10. **`MHRF-CORE-LSS-01`: LSS Fixed Regressor Precomputation (Core)**
    *   **Task:** Create function `prepare_lss_fixed_components_core(A_lss_fixed_matrix, intercept_column_in_Alss_idx = NULL, lambda_ridge_Alss = 0.1)`
        *   Input: `A_lss_fixed_matrix` (`n x q_lss` R matrix), `intercept_column_in_Alss_idx` (optional int), `lambda_ridge_Alss` (scalar).
        *   Output: List containing `P_lss_matrix` (`q_lss x n`) and `p_lss_vector` (`n x 1`).
        *   Details: Compute `P_lss` with ridge and `p_lss_vector` from `MASS::ginv(A_lss_fixed_matrix)`. Add jitter for matrix inversion.
    *   **DoD:** `P_lss_matrix` and `p_lss_vector` returned.

11. **`MHRF-CORE-LSS-02`: Voxel-Specific HRF Shape Reconstruction (Core)**
    *   **Task:** Create function `reconstruct_hrf_shapes_core(B_reconstructor_matrix, Xi_smoothed_matrix)`
        *   Input: `B_reconstructor_matrix` (`p x m`), `Xi_smoothed_matrix` (`m x V`).
        *   Output: `H_shapes_allvox_matrix` (`p x V` R matrix).
    *   **DoD:** Matrix of voxel-specific HRF shapes computed.

12. **`MHRF-CORE-LSS-03`: Implement Core Woodbury LSS Kernel (Per Voxel) (Core)**
    *   **Task:** Create internal function `run_lss_for_voxel_core(Y_proj_voxel_vector, X_trial_onset_list_of_matrices, H_shape_voxel_vector, A_lss_fixed_matrix, P_lss_matrix, p_lss_vector)`
        *   Input: `Y_proj_voxel_vector` (`n x 1` R vector), `X_trial_onset_list_of_matrices` (list of `T` `n x p` R matrices for trial onsets), `H_shape_voxel_vector` (`p x 1` R vector for HRF shape), `A_lss_fixed_matrix`, `P_lss_matrix`, `p_lss_vector`.
        *   Output: `beta_trial_vector` (`T x 1` R vector for one voxel).
        *   Details: Forms `C_vx` matrix, then Woodbury steps for `S_vx`, then `beta_trial_vx`.
    *   **DoD:** Function returns `T x 1` trial betas for one voxel. Tested.

13. **`MHRF-CORE-LSS-04`: Main LSS Loop and Optional `R_t` Precomputation (Core)**
    *   **Task:** Create function `run_lss_full_core(Y_proj_matrix, X_trial_onset_list_of_matrices, H_shapes_allvox_matrix, A_lss_fixed_matrix, P_lss_matrix, p_lss_vector, ram_heuristic_GB = 2)`
        *   Input: `Y_proj_matrix` (`n x V`), `X_trial_onset_list_of_matrices`, `H_shapes_allvox_matrix`, `A_lss_fixed_matrix`, `P_lss_matrix`, `p_lss_vector`, `ram_heuristic_GB` (for `R_t` precomputation).
        *   Output: `Beta_trial_allvox_matrix` (`T x V` R matrix).
        *   Details: Optionally precompute all `R_t_allvox` matrices. Loop over voxels (columns of `Y_proj_matrix`), calling `run_lss_for_voxel_core`.
    *   **DoD:** `Beta_trial_allvox_matrix` computed.

**EPIC 4: Core Engine - Alternating Optimization & Final Betas (Component 4)**

14. **`MHRF-CORE-ALTOPT-01`: Final Condition Beta Re-estimation (Core)**
    *   **Task:** Create function `estimate_final_condition_betas_core(Y_proj_matrix, X_condition_list_proj_matrices, H_shapes_allvox_matrix, lambda_beta_final, control_alt_list = list(max_iter=1, rel_change_tol=1e-4))`
        *   Input: `Y_proj_matrix` (`n x V`), `X_condition_list_proj_matrices` (list of `k` `n x p`), `H_shapes_allvox_matrix` (`p x V`), `lambda_beta_final` (scalar), `control_alt_list`.
        *   Output: `Beta_condition_final_matrix` (`k x V` R matrix).
        *   Details: Voxel-wise GLM. For each voxel `v`, construct design `X_design_v_final` by convolving each `X_condition_list_proj_matrices[[c]]` with `H_shapes_allvox_matrix[,v]`. Solve for betas. (MVP: `max_iter=1`).
    *   **DoD:** `Beta_condition_final_matrix` computed.

**EPIC 5: Neuroimaging Layer - Wrappers & I/O (MVP)**

15. **`MHRF-NIM-IO-MANIFOLD-01`: Manifold Construction Wrapper (Neuroimaging Layer)**
    *   **Task:** `construct_hrf_manifold_nim(hrf_library_source, ...)`
        *   Input: `hrf_library_source` (e.g., path to file with HRFs, list of `fmrireg::HRF` objects to be sampled into a matrix).
        *   Output: List with `B_reconstructor_matrix` and other manifold params.
        *   Details: Loads/generates `L_library_matrix`. Calls `MHRF-CORE-MANIFOLD-01` & `-02`.
    *   **DoD:** Wrapper produces manifold components from neuroimaging-relevant HRF sources.

16. **`MHRF-NIM-WRAP-SUBJECT-01`: Subject-Level Processing Wrapper (Neuroimaging Layer)**
    *   **Task:** `process_subject_mhrf_lss_nim(bold_input, mask_input, event_input, ...)`
        *   Input: `bold_input` (e.g., NIfTI path or `neuroim2::NeuroVec`), `mask_input` (path or `neuroim2::LogicalNeuroVol`), `event_input` (path to event file or `data.frame`), `manifold_objects`, `params_list`.
        *   Output: A list of R matrices (e.g., `Xi_smoothed_matrix`, `Beta_trial_matrix`).
        *   Details:
            *   Uses `neuroim2` to load BOLD/mask. Extracts `Y_data_matrix`.
            *   Uses `fmrireg` (`fmrireg::sampling_frame`, `fmrireg::event_model`) to prepare design matrix components (`X_condition_list_of_matrices`, `X_trial_onset_list_of_matrices`). This involves:
                *   Creating `fmrireg::event_spec` for conditions using `p` delta functions as the "HRF". Then `fmrireg::design_matrix(event_spec)` gives the `n x (kp)` matrix, which can be split into `k` `n x p` matrices for `X_condition_list`.
                *   Similar process for `X_trial_onset_list`.
            *   Calls the core engine functions (`MHRF-CORE-...`) in sequence.
    *   **DoD:** Wrapper function orchestrates calls to core engine, returning matrix outputs.

17. **`MHRF-NIM-OUTPUT-01`: Results Packaging (Neuroimaging Layer)**
    *   **Task:** Create function `package_mhrf_results_nim(core_results_list, reference_space, mask_vol)`
        *   Input: `core_results_list` (matrices from `MHRF-NIM-WRAP-SUBJECT-01`), `reference_space` (`neuroim2::NeuroSpace` from mask/BOLD), `mask_vol` (`neuroim2::LogicalNeuroVol`).
        *   Output: An S3/S4 object containing results as `neuroim2::NeuroVol` or `neuroim2::NeuroVec` objects.
        *   Details: Converts output matrices (e.g., `Xi_smoothed_matrix`) into `NeuroVol` objects using the mask and reference space.
    *   **DoD:** Structured output object with `neuroim2` data available. Basic `print`, `summary`, `plot` methods for this object.

**EPIC 6: Validation, QC, and Documentation (MVP)**

18. **`MHRF-VALIDATE-SIM-01`: Synthetic Data Validation Script (Core)**
    *   **Task:** Create `run_mhrf_lss_simulation_core()` script.
        *   Generates synthetic `Y_matrix`, true HRF parameters (as matrices), true trial/condition betas (as matrices). Can use `fmrireg::simulate_bold_signal()` to generate data for event regressors then combine with synthetic confound data.
        *   Calls the *core* M-HRF-LSS pipeline functions with these matrices.
        *   Computes RMSE, correlation for outputs against ground truth matrices.
    *   **DoD:** Script runs, produces validation metrics.

19. **`MHRF-QC-REPORT-01`: Basic QC Report Generation (Neuroimaging Layer)**
    *   **Task:** R Markdown template for HTML report.
        *   Inputs: Output object from `MHRF-NIM-OUTPUT-01`.
        *   Plots: Manifold eigenvalues, reconstruction error. Example HRF shapes. Map of one `Xi_smoothed` component (using `neuroim2::plot`).
    *   **DoD:** Template generates HTML report with specified plots.

20. **`MHRF-TEST-UNIT-CORE-01`: Unit Tests for Core Functions**
    *   **Task:** `testthat` tests for key core engine functions.
    *   Focus: Input/output dimensions, algorithmic correctness on small matrices.
    *   **DoD:** Core functions are unit tested.

21. **`MHRF-DOC-MVP-01`: MVP Documentation & Example Script**
    *   **Task:** Document main neuroimaging wrapper functions and the output object. Create a simple example script demonstrating the full neuroimaging layer pipeline.
    *   **DoD:** Basic user docs and a runnable example available.

---

This plan emphasizes building the computational heart first, making it easier to test and verify the algorithms independently before integrating the neuroimaging-specific data handling.