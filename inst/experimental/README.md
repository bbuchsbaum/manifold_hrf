# Experimental Functions

This directory contains experimental functions that are not currently used in the manifoldhrf package but may be useful for future development.

## Taylor Series Utilities

The files `taylor_utils.R` and `test-taylor-series-utils.R` contain functions for Taylor series expansion and numerical derivatives:

- `numeric_derivative()` - Compute numerical derivatives using central differences
- `numeric_taylor_coefficients()` - Generate Taylor series coefficients numerically  
- `evaluate_taylor()` - Evaluate Taylor polynomials at given points

These functions were added in PR #81 "for testing" but are not currently integrated into the M-HRF-LSS pipeline.

### Potential Future Uses

1. **Smooth HRF interpolation** - Taylor series could provide smooth interpolation between HRF grid points in the manifold
2. **Derivative-based optimization** - Numerical derivatives could be used for gradient-based HRF optimization
3. **Local HRF approximation** - Taylor expansions could approximate HRF shapes locally in the manifold space
4. **Sensitivity analysis** - Derivatives could help analyze how HRF parameters affect the BOLD signal

### Why Moved to Experimental

- Not currently used in any pipeline components
- Performance concerns due to extensive validation in the current implementation
- Unclear integration path with the manifold-based approach

The functions are preserved here for potential future development or research purposes.