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

### User Interfaces

The package provides multiple interfaces for different use cases:

1. **Simple Interface** (`R/simple_interface.R`):
   - `mhrf_lss()` - Streamlined function for users with pre-computed design matrices
   - `mhrf_lss_parameters()` - Parameter helper with presets

2. **fmrireg Integration** (`R/mhrf_lss_interface.R`):
   - `fit_mhrf_lss()` - Main interface following fmrireg patterns
   - Accepts formula specifications and fmri_dataset objects
   - Supports various HRF library sources and processing strategies

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
- `run_lss_for_voxel_corrected()`: Reference implementation for LSS (replaces former Woodbury variant)
- `run_lss_voxel_loop_core()`: Main loop with optional RAM-based precomputation
  controlled by `ram_heuristic_GB_for_Rt`. Trial regressors are precomputed
  when `T * V * 8 / 1e9` is below this limit.

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
   - RAM heuristics for trial-wise precomputation via
     `ram_heuristic_GB_for_Rt` in `run_lss_voxel_loop_core`
   - Chunked processing support via `process_in_chunks()`

3. **Numerical Stability**:
   - Ridge regularization throughout
   - QR decomposition for confound projection
   - SVD with near-zero singular value handling
   - Automatic regularization in robust variants

4. **Robustness & Fallbacks**:
   - Graceful degradation cascade: Manifold → PCA → Canonical HRF
   - Outlier detection and downweighting
   - Adaptive spatial smoothing based on SNR
   - Parameter presets (conservative/balanced/aggressive)

### Class System

- **`mhrf_lss_result`**: Main result object with S3 methods
  - `print()`, `summary()`, `plot()` methods
  - Access to core matrices and neuroimaging objects
  - QC metrics and diagnostics

- **`manifold_hrf_fit`**: R6 class for manifold HRF fitting
  - Encapsulates manifold construction and HRF estimation
  - Methods for prediction and evaluation

## Implementation Status

### Completed
- **Core Engine**: All components (0-4) fully implemented ✓
- **Neuroimaging Wrappers**: Basic integration complete ✓
- **User Interfaces**: Simple and fmrireg-style interfaces ✓
- **Validation Framework**: Comprehensive simulation tools ✓
- **QC & Reporting**: HTML report generation with diagnostics ✓
- **Robustness Suite**: Fallbacks, outlier handling, adaptive methods ✓
- **Testing**: 722+ tests across all components ✓

### Recent Additions
- Corrected LSS implementation (`run_lss_for_voxel_corrected()`)
- Unified user-friendly interface via `fit_mhrf_lss()`
- Enhanced parameter management with presets
- Improved memory strategies for large datasets
- Parallel processing utilities

### TODO
- Full neuroim2 integration for neuroimaging wrappers
- Advanced visualization methods
- GPU acceleration for large-scale processing
- Additional HRF library presets

## Testing Strategy

- Unit tests for each core function with synthetic data
- Integration tests for full pipeline components
- Validation against known ground truth in simulations
- Tests include dimension checking, numerical accuracy, and edge cases
- Memory and performance benchmarks for optimization

### Test Coverage
- **Total tests**: 722+ tests implemented
- **Component coverage**:
  - Manifold construction: 27 tests
  - Voxel-wise fit: 325 tests (most comprehensive)
  - Spatial smoothing: 51 tests
  - Trial-wise LSS: 45+ tests
  - Alternating optimization: 41 tests
  - Neuroimaging wrappers: 44 tests
  - Validation simulation: 93 tests
  - QC report generation: 34 tests
  - Soundness improvements: 62 tests

## Common Development Tasks

### Adding a New Core Function
1. Place in appropriate `R/core_*.R` file
2. Use matrix inputs/outputs only
3. Add comprehensive input validation
4. Include unit tests in `tests/testthat/test-*.R`
5. Document with roxygen2

### Modifying Pipeline Components
1. Check dependencies in other components
2. Update both core and wrapper layers if needed
3. Ensure backward compatibility
4. Run full test suite before committing

### Debugging Tips
- Use `manifold_hrf_logger()` for structured logging
- Check matrix dimensions with validation helpers
- Use `browser()` strategically in core loops
- Leverage diagnostic plots in QC reports