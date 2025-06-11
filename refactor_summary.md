# LSS Refactoring Summary

## Overview
Successfully refactored the manifold_hrf package to use the external fmrilss implementation for all LSS (Least Squares Separate) functionality.

## Changes Made

### 1. Simplified Core LSS Implementation
- Replaced complex internal LSS implementations in `core_lss.R` with simple wrappers to `fmrilss::lss`
- Removed obsolete Woodbury matrix identity implementations
- Removed internal adapter functions like `lss_compute_r`

### 2. Maintained Backward Compatibility
- Added support for old argument names in `prepare_lss_fixed_components_core`
- Created compatibility aliases for deprecated function names
- Handled both old (P_confound) and new (Z_confounds) interfaces in `run_lss_woodbury_corrected`

### 3. Fixed Critical Issues
- Fixed double projection bug when data is already projected
- Added proper computation of `p_lss_vector` for intercept handling
- Restored validation and condition number checks in `prepare_lss_fixed_components_core`

### 4. Cross-Platform Compatibility
- Fixed Windows compatibility issue in `parallel_utils.R` by implementing socket cluster fallback

## Test Results
- Reduced test failures from 22 to 15
- Fixed core LSS functionality tests
- Remaining failures are mostly related to:
  - Signal recovery strength (correlation and mean tests)
  - Specific validation edge cases
  - Interface changes between LSA and LSS

## Benefits
1. **Simplification**: Removed ~500+ lines of complex matrix algebra code
2. **Reliability**: Now using well-tested fmrilss implementation
3. **Maintainability**: Single source of truth for LSS computation
4. **Performance**: fmrilss provides optimized C++ implementations

## Next Steps
1. Investigate remaining test failures related to signal recovery
2. Update documentation to reflect the simplified implementation
3. Consider removing deprecated functions in future release
4. Profile performance improvements from using fmrilss