# manifold_hrf

Manifold Methods for Hemodynamic Response Function Analysis

## Overview

This package provides manifold learning methods for analyzing hemodynamic response functions (HRF) in neuroimaging data. It implements dimensionality reduction techniques and manifold-based approaches for understanding HRF patterns and dynamics. The package includes a complete pipeline for manifold-guided HRF estimation and trial-wise deconvolution (M-HRF-LSS) with support for neuroimaging data formats.

## Key Features

- **Manifold-guided HRF estimation** using diffusion maps
- **Trial-wise deconvolution** with Least Squares Separate (LSS)
- The function `run_lss_for_voxel_corrected()` now provides the
  reference implementation replacing the former Woodbury variant
  `run_lss_for_voxel_core()`.
- **Spatial smoothing** of manifold coordinates
- **Comprehensive QC reporting** with diagnostic plots
- **Integration** with `neuroim2` and `fmrireg` packages
- **Modular design** with separate core computational engine and neuroimaging layer

## Installation

```r
# Install development version from GitHub
# devtools::install_github("bbuchsbaum/manifold_hrf")
```

## Development Status

This package is currently under active development. The implementation follows the sprint plan outlined in `data-raw/Manifold_hrf_sprint.md`.

### Current Implementation Status

- ✅ Package structure and documentation framework
- 🚧 Core computational engine (Components 0-4)
- 🚧 Neuroimaging layer wrappers  
- 🚧 QC reporting framework
- 🚧 Validation and testing suite

## Quick Start

```r
library(manifold_hrf)

# Basic manifold construction functions are implemented. Neuroimaging
# wrappers currently provide minimal functionality.
# After running the pipeline you can visualize a QC summary:
# plot_qc_summary(result)
```

## Diagnostics and QC

After running the pipeline, inspect the resulting object with
`summary(result)` to view QC metrics such as mean R² and HRF statistics.
Visual diagnostics can be produced with:

```r
plot(result, type = "diagnostic")
```

For a comprehensive HTML report of all metrics and plots, call
`generate_qc_report()`.

## Development Roadmap

See the detailed sprint plan in `data-raw/Manifold_hrf_sprint.md` for the complete development roadmap, including:

- **EPIC 0-4**: Core computational engine implementation
- **EPIC 5**: Neuroimaging layer wrappers
- **EPIC 6**: Validation, QC, and documentation

## Contributing

This package follows the development sprint outlined in the project documentation. Please refer to the sprint plan for current priorities and implementation guidelines.

## License

MIT License 