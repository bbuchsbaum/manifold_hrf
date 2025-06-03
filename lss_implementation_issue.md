# LSS Implementation Discrepancies in manifoldhrf

## Summary

The package contains multiple implementations of LSS (Least Squares Separate) that produce different results. Analysis reveals that some implementations are not actually performing trial-wise LSS but rather simultaneous estimation, and the Woodbury-based optimization doesn't correctly implement the Woodbury matrix identity.

## Background

LSS (Least Squares Separate) is a method for estimating trial-specific amplitudes in fMRI analysis. For each trial, it estimates that trial's amplitude while treating all other trials as nuisance regressors. This is different from simultaneous estimation where all trials are estimated together.

## Current Implementations

### 1. `run_lss_for_voxel_corrected_full` 
- **Location**: `R/core_lss.R:96-132`
- **Issue**: Actually performs simultaneous estimation, NOT trial-wise LSS
- **Evidence**: Test in `test-lss-correction-validation.R:126` expects it to match simultaneous estimation

### 2. `run_lss_woodbury_corrected`
- **Location**: `R/core_lss.R:147-210`
- **Issue**: Uses correct trial-wise approach but implementation details differ from ground truth
- **Evidence**: Test failures in `test-lss-separate.R` show ~25% average differences from manual implementation

### 3. `.compute_lss_betas` (helper function)
- **Location**: `R/core_lss.R:69-80`
- **Issue**: Implements a scaling-based approach that doesn't match the Woodbury identity
- **Problem**: Uses `alpha = (1 - p'c) / ||v||²` scaling instead of proper matrix inversion update

## Test Failures

### 1. `test-lss-correction-validation.R`
```
Failure: Corrected implementation should match simultaneous estimation
Difference: 0.488 (expected < 1e-04)
```
This test actually confirms that the "corrected" implementation is doing simultaneous estimation.

### 2. `test-lss-separate.R`
```
Failure: `beta_wood` not equal to `beta_manual`
4/4 mismatches (average diff: 0.0822)
```
The Woodbury implementation produces systematically different results.

### 3. `test-lss-loop-core.R`
```
Failure: `Beta_core` not equal to `Beta_manual`
15/15 mismatches (average diff: 0.277)
```
Large discrepancies when using the loop implementation.

## Mathematical Analysis

### True LSS (for each trial t):
```
min_β ||y - X_t β_t - Σ(j≠t) X_j β_j - Z γ||² + λ||β||²
```

### What's being implemented:
1. **Simultaneous estimation**: `min_β ||y - Σ_j X_j β_j - Z γ||² + λ||β||²`
2. **Incorrect Woodbury**: Using scaling factors instead of: `(A'A + cc')^(-1) = (A'A)^(-1) - (A'A)^(-1)cc'(A'A)^(-1) / (1 + c'(A'A)^(-1)c)`

## Root Causes

1. **Conceptual confusion**: Some implementations solve for all trials simultaneously rather than separately
2. **Mathematical error**: The Woodbury identity is not correctly implemented
3. **Inconsistent projection**: Different methods handle confound removal differently
4. **Missing documentation**: No clear specification of what each function should compute

## Proposed Solutions

### Option 1: Fix implementations to match LSS specification
- Implement correct Woodbury identity for `run_lss_woodbury_corrected`
- Rewrite `run_lss_for_voxel_corrected_full` to do actual trial-wise estimation
- Ensure consistent confound handling across all methods

### Option 2: Clarify intentions and update tests
- If simultaneous estimation is actually desired for some methods, rename functions appropriately
- Update documentation to clearly state what each function computes
- Adjust tests to match actual intentions

### Option 3: Provide both approaches
- Clearly separate functions for:
  - `run_lss_simultaneous()` - current behavior of corrected_full
  - `run_lss_separate()` - true trial-wise LSS
  - `run_lss_woodbury()` - optimized version with correct identity

## Recommendations

1. **Immediate**: Document current behavior clearly in function headers
2. **Short-term**: Decide whether to fix implementations or adjust expectations
3. **Long-term**: Implement comprehensive tests comparing all methods on known ground truth

## References

- Original LSS paper: Turner et al. (2012) NeuroImage
- Woodbury matrix identity: https://en.wikipedia.org/wiki/Woodbury_matrix_identity
- Related issue in fmrireg: Consider how this integrates with the broader framework

## Code Examples

### Current problematic implementation:
```r
# From .compute_lss_betas - this doesn't implement LSS correctly
alpha_v_row <- (1 - pc_v_row) / pmax(cv_v_row, .Machine$double.eps)
S_effective_regressors_v <- sweep(V_regressors_v, MARGIN = 2,
                                  STATS = alpha_v_row, FUN = "*")
```

### What it should be (conceptually):
```r
# For each trial separately
for (t in 1:n_trials) {
  X_full <- cbind(X_trial[[t]], X_other_trials, Z_confounds)
  beta_all <- solve(crossprod(X_full) + lambda*I, crossprod(X_full, y))
  beta_trial[t] <- beta_all[1]  # Extract trial-specific estimate
}
```