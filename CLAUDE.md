# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

The `manifoldhrf` R package implements the Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS) pipeline for fMRI analysis. It uses manifold learning to estimate voxel-specific hemodynamic response functions (HRFs) and provides accurate trial-wise activation estimates.

## Key Commands

### Building and Testing
```bash
# Install the package
R CMD INSTALL .

# Run all tests
R -e "devtools::test()"

# Run specific test file
R -e "testthat::test_file('tests/testthat/test-manifold-construction.R')"

# Generate documentation
R -e "roxygen2::roxygenise()"

# Check package
R CMD check .

# Build package
R CMD build .
```

### Development Workflow
```bash
# Load package for interactive development
R -e "devtools::load_all()"

# Run specific function tests
R -e "library(manifoldhrf); library(testthat); test_file('tests/testthat/test-voxelwise-fit.R')"

# Lint code (if configured)
R -e "lintr::lint_package()"
```

## Architecture Overview

### Core vs Neuroimaging Layer Separation

The package follows a strict separation between:

1. **Core computational engine** (`R/core_*.R`): 
   - Operates on standard R matrices and lists
   - No neuroimaging dependencies
   - Fully testable with synthetic data
   - Matrix/sparse matrix operations using base R, Matrix, and RSpectra packages

2. **Neuroimaging wrapper layer** (`R/neuroimaging_*.R`):
   - Handles neuroim2/fmrireg data structures
   - NIfTI I/O operations
   - Converts between neuroimaging formats and core matrix inputs

### M-HRF-LSS Pipeline Components

The pipeline consists of 5 main components (0-4) plus pre-flight QC (-1):

**Component 0: HRF Manifold Construction** (Run once per HRF library)
- `calculate_manifold_affinity_core()`: Creates Markov transition matrix using self-tuning local scaling
- `get_manifold_basis_reconstructor_core()`: Computes diffusion map embedding and reconstructor matrix

**Component 1: Voxel-wise HRF Estimation**
- `project_out_confounds_core()`: QR-based confound removal from data and designs
- `transform_designs_to_manifold_basis_core()`: Projects design matrices to manifold space
- `solve_glm_for_gamma_core()`: Ridge regression for gamma coefficients with optional orthogonal approximation
- `extract_xi_beta_raw_svd_core()`: SVD-based extraction of manifold coordinates and amplitudes
- `apply_intrinsic_identifiability_core()`: Sign/scale constraints for identifiability

**Component 2: Spatial Smoothing**
- `make_voxel_graph_laplacian_core()`: k-NN graph construction for spatial regularization
- `apply_spatial_smoothing_core()`: Graph Laplacian-based smoothing of manifold coordinates

**Component 3: Trial-wise LSS**
- `prepare_lss_fixed_components_core()`: Precompute fixed regression components
- `reconstruct_hrf_shapes_core()`: Transform smoothed manifold coords back to HRF shapes
- `run_lss_for_voxel_core()`: Woodbury-optimized LSS for single voxel
- `run_lss_voxel_loop_core()`: Main loop with optional RAM-based precomputation

**Component 4: Alternating Optimization**
- `estimate_final_condition_betas_core()`: Re-estimate condition betas using final HRFs

### Key Design Patterns

1. **Matrix Organization**: 
   - Data matrices: `n x V` (timepoints × voxels)
   - HRF library: `p x N` (HRF samples × number of HRFs)
   - Design matrices: `n x p` (timepoints × HRF length)
   - Manifold coordinates: `m x V` (dimensions × voxels)

2. **Memory Management**:
   - Sparse matrix support for large HRF libraries (N > 5000)
   - Optional HDF5 backing for very large datasets
   - RAM heuristics for trial-wise precomputation

3. **Numerical Stability**:
   - Ridge regularization throughout
   - QR decomposition for confound projection
   - SVD with near-zero singular value handling

## Implementation Status

### Completed
- Component 0: HRF Manifold Construction (2 functions) ✓
  - `calculate_manifold_affinity_core()` - Creates Markov transition matrix
  - `get_manifold_basis_reconstructor_core()` - Computes diffusion map and reconstructor
- Component 1: Voxel-wise Manifold Fit (5 functions) ✓
  - `project_out_confounds_core()` - QR-based confound removal
  - `transform_designs_to_manifold_basis_core()` - Transform designs to manifold space
  - `solve_glm_for_gamma_core()` - Ridge regression for gamma coefficients
  - `extract_xi_beta_raw_svd_core()` - SVD extraction of manifold coords and amplitudes
  - `apply_intrinsic_identifiability_core()` - Sign/scale constraints for identifiability
- Component 2: Spatial Smoothing (2 functions) ✓
  - `make_voxel_graph_laplacian_core()` - k-NN graph Laplacian construction
  - `apply_spatial_smoothing_core()` - Laplacian regularization for spatial smoothing

- Component 3: Trial-wise LSS (4 functions) ✓
  - `prepare_lss_fixed_components_core()` - Precompute fixed LSS components
  - `reconstruct_hrf_shapes_core()` - Transform manifold coords to HRF shapes
  - `run_lss_for_voxel_core()` - Woodbury LSS for single voxel
  - `run_lss_voxel_loop_core()` - Main LSS loop with RAM optimization
- Component 4: Alternating Optimization (1 function) ✓
  - `estimate_final_condition_betas_core()` - Re-estimate condition betas with final HRFs

- **Neuroimaging Wrappers** (partially implemented)
  - `construct_hrf_manifold_nim()` ✓ - Enhanced manifold construction with fmrireg integration
  - `process_subject_mhrf_lss_nim()` (placeholder) - Subject-level processing wrapper
  - `package_mhrf_results_nim()` (placeholder) - Results packaging with neuroim2 integration

- **Validation Framework** ✓
  - `run_mhrf_lss_simulation()` - Main simulation runner with noise robustness testing
  - Helper functions for generating ground truth HRFs, designs, and amplitudes
  - Comprehensive evaluation metrics for HRF recovery and amplitude estimation
  - Noise robustness curve generation

- **QC Report Generation** ✓
  - `generate_qc_report()` - Creates comprehensive HTML QC reports with R Markdown
  - `create_qc_flags()` - Evaluates pipeline results and generates QC flags
  - `extract_hrf_stats()` - Computes HRF shape statistics (peak time, FWHM, etc.)
  - `compute_qc_diagnostics()` - Calculates diagnostic metrics (R², reconstruction error)
  - Comprehensive report template with manifold diagnostics, HRF stats, and failure mode checklist

- **Soundness & Robustness Improvements** ✓
  - `get_manifold_basis_reconstructor_robust()` - Manifold construction with PCA fallback
  - `check_hrf_library_quality()` - Validates HRF library quality
  - `suggest_parameters()` - Data-driven parameter selection
  - `get_preset_params()` - Conservative/balanced/aggressive presets
  - `extract_xi_beta_raw_svd_robust()` - Robust SVD with automatic regularization
  - `apply_spatial_smoothing_adaptive()` - SNR-adaptive spatial smoothing
  - `detect_outlier_timepoints()` - Outlier detection and downweighting
  - `screen_voxels()` - Quality-based voxel screening
  - `run_with_fallback_cascade()` - Graceful degradation (Manifold → PCA → Canonical)
  - `estimate_memory_requirements()` - Memory usage estimation
  - `process_in_chunks()` - Memory-efficient chunked processing

### TODO
- Complete neuroimaging wrapper implementations (requires neuroim2/fmrireg)
- Pipeline orchestration functions
- Performance optimization (parallelization, C++ acceleration)

## Testing Strategy

- Unit tests for each core function with synthetic data
- Integration tests for full pipeline components
- Validation against known ground truth in simulations
- Tests include dimension checking, numerical accuracy, and edge cases

### Test Coverage
- **Total tests**: 722 tests implemented
- **Component coverage**:
  - Manifold construction: 27 tests
  - Voxel-wise fit: 325 tests (most comprehensive)
  - Spatial smoothing: 51 tests
  - Trial-wise LSS: 45 tests
  - Alternating optimization: 41 tests
  - Neuroimaging wrappers: 44 tests
  - Validation simulation: 93 tests
  - QC report generation: 34 tests
  - Soundness improvements: 62 tests