Okay, this is an excellent synthesis of the iterative discussions. I will now revise the original `Manifold.txt` proposal to incorporate all the refinements, decisions, and new components discussed. The goal is a single, comprehensive engineering blueprint reflecting the final state of our design consensus.

---
**Revised Proposal: Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS) v3.0**

**1. Executive Summary:**

This proposal outlines an advanced fMRI analysis pipeline, "Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS)." It combines a novel HRF estimation technique based on manifold learning with a highly optimized Least Squares Separate (LSS) kernel for single-trial beta estimation. The M-HRF-LSS pipeline aims to provide:
1.  Spatially smooth and physiologically plausible voxel-wise HRF estimates.
2.  Associated per-condition activation amplitudes, refined after HRF estimation.
3.  Accurate per-trial activation amplitudes.
4.  A comprehensive QC report and validation metrics for scientific credibility.
All components are designed for computational efficiency, primarily relying on closed-form solutions, BLAS operations, sparse matrix algebra, and optional C++ acceleration, making it suitable for large-scale whole-brain analyses on standard CPU hardware, with considerations for memory management via HDF5-backed arrays for very large datasets.

**2. Core Components & Algorithmic Workflow:**

**Component -1: Pre-flight Data Quality Control (QC)**

*   **Objective:** Perform initial checks on input data to ensure basic validity and flag potential issues.
*   **Inputs:**
    *   Raw BOLD data (e.g., NIfTI file path or `n x V` matrix).
    *   Trial/event timing files.
    *   Motion parameter files.
    *   TR (Repetition Time).
*   **Steps:**
    1.  **TR Uniformity:** Check if TR is consistent.
    2.  **Motion Spikes:** Identify volumes with excessive motion (e.g., >2mm framewise displacement, using tools like `fsl_motion_outliers`).
    3.  **Trial Count:** Minimum trial count per condition (e.g., warn if < 10).
    4.  **Trial Density:** Warn if trial density is very low (e.g., < 0.1 trials per TR * run length).
    5.  **Noise Level:** Estimate mean DVARS; warn if excessively high (e.g., > 5%).
*   **Outputs:** A list of `qc_flags` (e.g., `low_trial_count`, `high_motion_spikes`, `high_noise_dvars`) passed to subsequent components and included in the final report.

**Component 0: HRF Manifold Construction (Once per Study/HRF Library)**

*   **Objective:** Create a low-dimensional, smooth representation (manifold) of plausible HRF shapes.
*   **Inputs:**
    *   `L_library`: `p x N` matrix of `N` candidate HRF shapes, each sampled at `p` time points (e.g., from a physiological model grid, half-cosine set, FLOBS).
    *   `m_manifold_dim_target`: Target dimensionality of the manifold (user suggestion, e.g., 3-5).
    *   `m_manifold_dim_min_variance`: Minimum variance to be explained by manifold (e.g., 0.95 for automatic `m` selection).
    *   `k_local_nn_for_sigma`: k-Nearest Neighbors for self-tuning bandwidth `σ_i` (default=7; user-exposed).
    *   (Optional) `distance_engine`: Method for distance calculation (default="euclidean"; options: "euclidean", "ann_euclidean" via RcppHNSW if N > 10k, or custom plugin like "geodesic").
    *   (Optional) `sparse_if_N_gt`: Threshold for N to switch to sparse affinity matrix (default=5000).
    *   (Optional) `k_nn_for_W_sparse`: Number of neighbors for sparse affinity matrix `W` if sparsified (default=20).
*   **Steps:**
    1.  **Affinity Matrix `W` (Self-Tuning):**
        *   Compute pairwise distances `dists_ij` between all HRFs in `L_library` using `distance_engine`.
        *   For each HRF `i`, find its `k_local_nn_for_sigma`-th nearest neighbor distance `σ_i`.
        *   `W_ij = exp(-dists_ij² / (σ_i * σ_j))` (Local scaling heuristic, Zelnik-Manor & Perona, 2005).
        *   If `N > sparse_if_N_gt`, `W` can be sparsified (e.g., k-NN graph using `k_nn_for_W_sparse`).
    2.  **Markov Matrix `S`:** `D_inv = diag(1/rowSums(W))`, `S = D_inv %*% W`.
    3.  **Diffusion Map Coordinates `Φ_raw` & Automatic `m` selection:**
        *   Compute top `k_eig_max` (e.g., 10 or `m_manifold_dim_target + 5`) eigenvectors of `S` (e.g., via `RSpectra::eigs_sym`). `Φ_raw_full = eig_S$vectors`.
        *   `eigenvalues_S = eig_S$values`.
        *   `cum_var_explained = cumsum(eigenvalues_S[-1]) / sum(eigenvalues_S[-1])` (excluding trivial eigenvector).
        *   `m_auto = which(cum_var_explained >= m_manifold_dim_min_variance)[1]`.
        *   `m_manifold_dim = min(m_auto, m_manifold_dim_target)` if user specified target, else `m_auto`. Warn if user-set `m_manifold_dim` is much lower than `m_auto`.
        *   Record `m_auto` and chosen `m_manifold_dim`.
    4.  **Manifold Coordinates `Φ`:** `Φ = Φ_raw_full[, 2:(m_manifold_dim + 1)]` (an `N x m_manifold_dim` matrix). Enforce consistent sign (e.g., first element of each eigenvector positive) for reproducibility of `B_reconstructor`.
    5.  **HRF Reconstructor `B_reconstructor`:** `B_reconstructor = L_library %*% Φ %*% solve(crossprod(Φ) + 1e-8 * diag(m_manifold_dim))` (a `p x m_manifold_dim` matrix). This maps `m`-dimensional manifold coordinates `ξ` back to a `p`-sample HRF shape: `h_shape = B_reconstructor %*% ξ`.
*   **Outputs:** `B_reconstructor`, `Φ`, `eigenvalues_S`, `m_manifold_dim`, `m_auto`, parameters used.

**Component 1: Voxel-wise HRF Manifold Coordinate & Initial Condition Amplitude Estimation**

*   **Objective:** For each voxel, estimate its `m`-dimensional HRF manifold coordinate `ξ_v` and `k` per-condition amplitudes `β_v`.
*   **Inputs:**
    *   `Y_bold`: `n x V` BOLD data matrix.
    *   `X_condition_list`: List of `k` matrices, where `X_condition_list[[c]]` is an `n x p` Toeplitz design matrix for condition `c` (convolving onsets with `p` delta functions).
    *   `Z_confounds`: (Optional) `n x q_confound` matrix of nuisance regressors.
    *   `B_reconstructor`: `p x m` matrix from Component 0.
    *   `lambda_gamma`: Ridge penalty for the `(km)x(km)` GLM solve.
    *   `h_ref_shape_canonical`: `p x 1` canonical HRF shape for identifiability.
    *   `ident_scale_method`: Method for scaling reconstructed HRFs ("l2_norm", "max_abs_val", "none"). Default "l2_norm".
    *   `ident_sign_method`: Method for sign alignment ("canonical_correlation", "data_fit_correlation"). Default "canonical_correlation".
    *   `orthogonal_approx_flag`: Boolean (default `FALSE`).
*   **Steps:**
    1.  **Confound Projection (if `Z_confounds` provided):**
        *   `Q_Z = qr.Q(qr(Z_confounds, LAPACK=TRUE))`.
        *   `Y_proj = Y_bold - Q_Z %*% tcrossprod(Q_Z, Y_bold)`.
        *   `X_condition_list_proj = lapply(X_condition_list, function(X) X - Q_Z %*% tcrossprod(Q_Z, X))`.
    2.  **Per-Condition Design in Manifold Basis `Z_list`:**
        *   `Z_list = lapply(X_condition_list_proj, function(Xc) Xc %*% B_reconstructor)`. Each `Z_list[[c]]` is `n x m`.
    3.  **Combined Design `X_tilde` and GLM Solve for `Gamma_coeffs`:**
        *   `X_tilde = do.call(cbind, Z_list)` (an `n x (km)` matrix).
        *   `XtX_tilde_reg = crossprod(X_tilde) + lambda_gamma * diag(k*m)`. (Consider block-wise computation if `X_tilde` is too large for memory).
            *   If `orthogonal_approx_flag = TRUE`, zero out off-diagonal `m x m` blocks of `crossprod(X_tilde)` before adding ridge.
        *   `XtY_tilde = crossprod(X_tilde, Y_proj)`.
        *   `Gamma_coeffs = solve(XtX_tilde_reg, XtY_tilde)` (a `(km) x V` matrix).
    4.  **Extract `ξ_raw_allvox` and `β_raw_allvox` via SVD (per voxel):**
        *   Initialize `Xi_raw_allvox = matrix(0, m, V)` and `Beta_raw_allvox = matrix(0, k, V)`.
        *   Loop `vx` from 1 to `V`: (Consider power iteration if `m` is large and `k` small).
            *   `G_vx = matrix(Gamma_coeffs[, vx], nrow = m, ncol = k)`.
            *   `svd_G_vx = svd(G_vx)`. Handle near-zero singular values (set `xi_vx`, `beta_vx` to zero).
            *   `Xi_raw_allvox[,vx] = svd_G_vx$u[,1] * sqrt(svd_G_vx$d[1])`.
            *   `Beta_raw_allvox[,vx] = svd_G_vx$v[,1] * sqrt(svd_G_vx$d[1])`.
    5.  **Intrinsic Identifiability (on `ξ` and `β`):**
        *   `xi_ref_coord = MASS::ginv(B_reconstructor) %*% h_ref_shape_canonical`.
        *   For each voxel `vx`:
            *   `xi_vx = Xi_raw_allvox[,vx]`, `beta_vx = Beta_raw_allvox[,vx]`.
            *   `reconstructed_hrf_vx_unscaled = B_reconstructor %*% xi_vx`.
            *   **Sign:**
                *   If `ident_sign_method == "canonical_correlation"`: `sgn = sign(sum(xi_vx * xi_ref_coord))`.
                *   If `ident_sign_method == "data_fit_correlation"`: (Requires original `Y_proj` and `X_condition_list_proj` for this voxel) Estimate full model fit for `+sgn` and `-sgn` choice for `xi_vx`, choose sign that maximizes `R^2`.
            *   `xi_vx_signed = xi_vx * sgn`, `beta_vx_signed = beta_vx * sgn`.
            *   `reconstructed_hrf_vx_signed = B_reconstructor %*% xi_vx_signed`.
            *   **Scale:**
                *   If `ident_scale_method == "l2_norm"`: `scl = 1 / pmax(sqrt(sum(reconstructed_hrf_vx_signed^2)), .Machine$double.eps)`.
                *   If `ident_scale_method == "max_abs_val"`: `scl = 1 / pmax(max(abs(reconstructed_hrf_vx_signed)), .Machine$double.eps)`.
                *   Else (`"none"`): `scl = 1`.
            *   `Xi_ident_allvox[,vx] = xi_vx_signed * scl`.
            *   `Beta_ident_allvox[,vx] = beta_vx_signed / scl`. (Ensure zeroed if scl was based on machine epsilon).
*   **Outputs:** `Xi_ident_allvox` (`m x V`), `Beta_ident_allvox` (`k x V` initial condition betas).

**Component 2: Spatial Smoothing of Manifold Coordinates**

*   **Objective:** Regularize HRF manifold coordinates spatially across the brain.
*   **Inputs:**
    *   `Xi_ident_allvox`: `m x V` matrix from Component 1.
    *   `voxel_coordinates_matrix`: `V x 3` matrix for graph construction.
    *   `lambda_spatial_smooth`: Spatial smoothing strength.
    *   (Optional) `smoothing_engine`: Method for graph Laplacian and solve (default="knn_graph_laplacian"; options: "knn_graph_laplacian", "surface_geodesic_laplacian").
    *   (Optional) `num_neighbors_Lsp`: Number of neighbors for k-NN graph Laplacian (default=6; 6 for face-connected up to 26 for corner-connected).
*   **Steps:**
    1.  **Construct `voxel_graph_laplacian_Lsp` (`V x V` sparse matrix) using `voxel_coordinates_matrix`, `num_neighbors_Lsp` and `smoothing_engine`.** (Consider RcppAnnoy or similar for k-NN search if pure R is slow). MVP: volume-based k-NN.
    2.  Initialize `Xi_smoothed_allvox = matrix(0, m, V)`.
    3.  For each manifold dimension `j` from 1 to `m`:
        *   `Xi_smoothed_allvox[j,] = Matrix::solve(A = Diagonal(V) + lambda_spatial_smooth * voxel_graph_laplacian_Lsp, b = Xi_ident_allvox[j,])`.
*   **Outputs:** `Xi_smoothed_allvox` (`m x V`). Helper: L-curve or GCV for `lambda_spatial_smooth` selection.

**Component 3: Trial-wise Amplitude Estimation (LSS using Smoothed Manifold HRFs)**

*   **Objective:** Estimate single-trial amplitudes `β_trial` using the spatially smoothed, voxel-specific HRFs.
*   **Inputs:**
    *   `Y_proj`: `n x V` confound-projected BOLD data.
    *   `X_trial_onset_list`: List of `T` matrices, where `X_trial_onset_list[[t]]` is an `n x p` Toeplitz design for trial `t`'s onset. (Consider `upsample_factor` for sparse designs).
    *   `B_reconstructor`: `p x m` from Component 0.
    *   `Xi_smoothed_allvox`: `m x V` from Component 2.
    *   `A_lss_fixed`: `n x q_lss` matrix of fixed nuisance regressors for LSS (e.g., *original* `Z_confounds` if not task-related, or just an intercept).
    *   `intercept_column_in_Alss`: Index of intercept column in `A_lss_fixed` if present, else `NULL`.
    *   `lambda_ridge_Alss`: Ridge penalty for `(A_lss_fixedᵀA_lss_fixed + λI)⁻¹A_lss_fixedᵀ`.
*   **Steps (Voxel-Batched Woodbury - Mode B from LSS proposal):**
    1.  **Precompute LSS Fixed Components:**
        *   `P_lss_AtA = (crossprod(A_lss_fixed) + lambda_ridge_Alss * diag(ncol(A_lss_fixed)))`. Add tiny jitter (e.g., `1e-6 * median(diag(P_lss_AtA))`) to diagonal if near-singular. Warn user.
        *   `P_lss = solve(P_lss_AtA, t(A_lss_fixed))` (`q_lss x n`).
        *   `p_lss_vec`: If `intercept_column_in_Alss` is not `NULL`, `p_lss_vec = ginv(A_lss_fixed)[intercept_column_in_Alss,]`. Else `p_lss_vec = rep(0, n)`.
    2.  **Precompute Voxel-Specific HRF Shapes `H_shapes_allvox`:**
        *   `H_shapes_allvox = B_reconstructor %*% Xi_smoothed_allvox` (a `p x V` matrix).
    3.  **Optional (RAM permitting via heuristic: `T*V*8 < 0.5 * available_RAM`): Precompute all `R_t_allvox` matrices:**
        *   For trial `t = 1...T`: `R_t_allvox = X_trial_onset_list[[t]] %*% H_shapes_allvox` (an `n x V` matrix).
    4.  **Initialize `Beta_trial_allvox = matrix(0, T, V)`.**
    5.  **Voxel Loop (parallelizable, consider RcppArmadillo w/ OpenMP if slow):** For each voxel `v = 1...V`:
        a.  **Construct `C_v` (`n x T`):** `C_v[,t] = (if R_t precomputed) R_t_allvox[,v] else (X_trial_onset_list[[t]] %*% H_shapes_allvox[,v])`.
        b.  **Woodbury LSS Core (for current voxel `v`):**
            i.  `U_v = P_lss %*% C_v` (`q_lss x T`).
            ii. `V_regressors_v = C_v - A_lss_fixed %*% U_v` (`n x T`).
            iii. `pc_v_row = crossprod(p_lss_vec, C_v)` (`1 x T`).
            iv. `cv_v_row = colSums(V_regressors_v * V_regressors_v)` (`1 x T`).
            v.  `alpha_v_row = (1 - pc_v_row) / pmax(cv_v_row, .Machine$double.eps)`.
            vi. `S_effective_regressors_v = sweep(V_regressors_v, MARGIN = 2, STATS = alpha_v_row, FUN = "*")`.
            vii. `S_effective_regressors_v = sweep(S_effective_regressors_v, MARGIN = 1, STATS = p_lss_vec, FUN = "+")`.
        c.  **Betas for Voxel `v`: `Beta_trial_allvox[,v] = crossprod(S_effective_regressors_v, Y_proj[,v])` (`T x 1`).**
*   **Outputs:** `Beta_trial_allvox` (`T x V`).

**Component 4: Alternating Optimization & Final Condition Beta Re-estimation (Optional)**

*   **Objective:** Refine condition-level betas using the spatially smoothed HRFs.
*   **Inputs:**
    *   `Y_proj`: `n x V` confound-projected BOLD data.
    *   `X_condition_list_proj`: List of `k` projected design matrices from Component 1.
    *   `H_shapes_allvox`: `p x V` smoothed HRF shapes from Component 3.
    *   `lambda_beta_final`: Ridge penalty for final beta estimation.
    *   `control_alt`: List with `max_iter` (e.g., 5) and `rel_change_tol` (e.g., 1e-4).
*   **Steps (Iterative refinement, default 1 pass for beta re-estimation):**
    1.  Initialize `Beta_condition_final_allvox = matrix(0, k, V)`.
    2.  **Loop `iter` from 1 to `control_alt$max_iter` (default `max_iter=1` for MVP, effectively just re-fitting betas once):**
        a.  **Construct Condition-Specific Design with Smoothed HRFs per voxel and condition:**
            For each voxel `v` and condition `c`, form `X_design_v_final[,c] = X_condition_list_proj[[c]] %*% H_shapes_allvox[,v]`.
        b.  **Voxel Loop:** For each voxel `v`:
            i.  `X_design_v_final = matrix(NA, nrow=n, ncol=k)`.
            ii. For condition `c = 1...k`: `X_design_v_final[,c] = X_condition_list_proj[[c]] %*% H_shapes_allvox[,v]`.
            iii. `XtX_v_final_reg = crossprod(X_design_v_final) + lambda_beta_final * diag(k)`.
            iv. `XtY_v_final = crossprod(X_design_v_final, Y_proj[,v])`.
            v.  `Beta_condition_final_allvox[,v] = solve(XtX_v_final_reg, XtY_v_final)`.
        c.  (If `max_iter > 1`): Re-estimate `Xi_smoothed_allvox` using current `Beta_condition_final_allvox` (similar to Component 1 but fixing betas and solving for Xi, then re-smooth). Check for convergence (`rel_change(Beta_condition_final_allvox)`).
*   **Outputs:** `Beta_condition_final_allvox` (`k x V`).

**3. Final Outputs of M-HRF-LSS Pipeline (Stored in an S3/S4 Object):**

*   `Parameters`: All input parameters, software versions, `qc_flags` from Component -1.
*   `Manifold`:
    *   `B_reconstructor`: `p x m` HRF reconstructor matrix.
    *   `Phi_library_coords`: `N x m` manifold coordinates of library HRFs.
    *   `eigenvalues_S`: Eigenvalues from diffusion map.
    *   `m_manifold_dim`, `m_auto_selected`.
*   `HRF_Estimates`:
    *   `Xi_ident_allvox`: `m x V` identifiability-constrained manifold coordinates (before spatial smoothing).
    *   `H_reconstructed_raw_allvox`: `p x V` HRFs from `Xi_ident_allvox`.
    *   `Xi_smoothed_allvox`: `m x V` spatially smoothed manifold coordinates.
    *   `H_reconstructed_smoothed_allvox`: `p x V` final reconstructed HRF shapes from `Xi_smoothed_allvox`.
*   `Amplitude_Estimates`:
    *   `Beta_condition_initial_allvox`: `k x V` condition-level amplitudes from Component 1.
    *   `Beta_condition_final_allvox`: `k x V` condition-level amplitudes from Component 4.
    *   `Beta_trial_allvox`: `T x V` trial-level LSS amplitudes.
*   `QC_Report_Path`: Path to the generated HTML QC report.
*   `qc_struct`: Programmatic access to QC flags (e.g. `qc_flags$low_trial_count`).

**4. Validation Framework (VAL-SIM-01 Integration)**

*   **Objective:** Provide tools to validate pipeline performance using synthetic data.
*   **Implementation:** A separate function/script `run_mhrf_lss_simulation()`:
    *   Generates synthetic BOLD data with known:
        *   HRF shapes (canonical, delayed, dispersed, biphasic).
        *   Activation locations and trial amplitudes.
        *   Noise levels (e.g., 0-10% DVARS, via additive white Gaussian noise).
    *   Runs the M-HRF-LSS pipeline on this synthetic data.
    *   Computes and reports metrics:
        *   **Temporal:** Peak-time error (s), FWHM error (s) for HRFs.
        *   **Spatial:** Dice coefficient on supra-threshold activation masks.
        *   **Amplitude:** Pearson ρ between true and estimated trial/condition betas.
        *   **Noise Robustness:** Performance-vs-noise curves for key metrics.
*   **Output:** A simulation report, figures, and summary statistics.

**5. Pipeline Quality Control (QC) Report (QC-REPORT-01 Integration)**

*   **Objective:** Generate a per-run HTML report summarizing key processing steps and diagnostics.
*   **Implementation:** An R Markdown template knitted at the end of the pipeline.
*   **Contents:**
    *   Input parameters and `qc_flags` from Component -1.
    *   **Manifold Diagnostics (DIAG-MANIFOLD-01):**
        *   Diffusion map eigenvalue spectrum (scree plot).
        *   Cumulative variance explained by `m` dimensions, chosen `m` vs `m_auto`.
        *   Reconstruction error (`||L - BΦ||F / ||L||F`) vs `m`.
        *   (Optional) Interactive 3D scatter plot of `Φ` if `m=3`.
    *   **HRF & Beta Diagnostics:**
        *   Maps of `Xi_smoothed_allvox` (e.g., first 3 components).
        *   Before/after spatial smoothing maps of one `ξ` coordinate.
        *   Histogram of HRF peak times and FWHM from `H_reconstructed_smoothed_allvox`.
        *   Median voxel-wise R² of GLM fit (Component 1) and LSS fit (Component 3).
        *   Convergence plot for alternating optimization (ΔHRF/ΔBeta per iteration) if Component 4 `max_iter > 1`.
    *   **QC Flags Summary:** `low_trial_count`, `poor_fits` (e.g., mean `R2_vox < 0.1` for >30% voxels), `unstable_hrf` (e.g., HRF peaks < 2s or > 10s), `high_motion_spikes`, `high_noise_dvars`. Each with red/yellow/green badge.
    *   **Failure-Mode Checklist:** Guidance on interpreting flags.
    *   Elapsed time and iteration counts for key components.

