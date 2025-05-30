# manifold_hrf

Manifold Methods for Hemodynamic Response Function Analysis

## Overview

This package provides manifold learning methods for analyzing hemodynamic response functions (HRF) in neuroimaging data. It implements dimensionality reduction techniques and manifold-based approaches for understanding HRF patterns and dynamics. The package includes a complete pipeline for manifold-guided HRF estimation and trial-wise deconvolution (M-HRF-LSS) with support for neuroimaging data formats.

## Key Features

- **Manifold-guided HRF estimation** using diffusion maps
- **Trial-wise deconvolution** with Least Squares Separate (LSS) 
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

This package is currently under active development. The implementation follows the sprint plan outlined in `data-raw/Mainfold_hrf_sprint.md`.

### Current Implementation Status

- ✅ Package structure and documentation framework
- 🚧 Core computational engine (Components 0-4)
- 🚧 Neuroimaging layer wrappers  
- 🚧 QC reporting framework
- 🚧 Validation and testing suite

## Quick Start

```r
library(manifold_hrf)

# This package is not yet functional - functions are placeholder stubs
# Please refer to the sprint plan for implementation roadmap
```

## Development Roadmap

See the detailed sprint plan in `data-raw/Mainfold_hrf_sprint.md` for the complete development roadmap, including:

- **EPIC 0-4**: Core computational engine implementation
- **EPIC 5**: Neuroimaging layer wrappers
- **EPIC 6**: Validation, QC, and documentation

## Contributing

This package follows the development sprint outlined in the project documentation. Please refer to the sprint plan for current priorities and implementation guidelines.

## License

MIT License 