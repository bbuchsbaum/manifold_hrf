# Test Failure Summary for manifoldhrf Package

## Overall Statistics
- **Total Tests**: 1,412 (1,275 passed + 121 warnings + 9 failures + 7 skipped)
- **Failed Tests**: 9
- **Test Files with Failures**: 4

## Detailed Failure Analysis

### 1. mhrf_lss pipeline integration (`test-mhrf-lss-pipeline.R`)
**Failures**: 2

#### Test: "LSS failure propagation does not crash pipeline" (Line 41)
- **Error**: `chol.default(XtX): the leading minor of order 1 is not positive`
- **Location**: `fmrilss::lss()` → `chol2inv(chol(XtX))`
- **Root Cause**: The matrix XtX is not positive definite, likely due to rank deficiency in the design matrix
- **Context**: Testing error handling when LSS fails

#### Test: "metadata alignment handles non sequential trials" (Line 141)
- **Error**: Same as above - `chol.default(XtX): the leading minor of order 1 is not positive`
- **Root Cause**: Similar issue with matrix singularity in the LSS computation

### 2. Core voxel fit engine stability (`test-voxelfit-engine-collinearity.R`)
**Failures**: 1

#### Test: "solve_glm_for_gamma_core produces unstable betas for near collinear design" (Line 34)
- **Error**: `qr.solve(XtX_reg, XtY): singular matrix 'a' in solve`
- **Location**: `solve_glm_for_gamma_core()` at line 130 of voxel_fit_core.R
- **Root Cause**: Near-collinear design matrix with lambda_gamma=0 causing numerical instability
- **Context**: Testing behavior with near-collinear designs

### 3. Voxelfit engine improved (`test-voxelfit-engine-improved.R`)
**Failures**: 1

#### Test: "apply_identifiability_vectorized is efficient" (Line 196)
- **Error**: `all(abs(beta_norms - 1) < 1e-10 | beta_norms == 0) is not TRUE`
- **Root Cause**: Beta normalization check failing - some beta values are not properly normalized to unit norm
- **Context**: Testing the vectorized identifiability constraints

### 4. Voxelwise fit (`test-voxelwise-fit.R`)
**Failures**: 5

#### Test: "extract_xi_beta_raw_svd_core validates inputs" (Lines 490, 495)
- **Error 1**: Expected "positive integer" but got "Gamma_coeffs_matrix has 6 rows but expected 0 (k * m)"
- **Error 2**: Expected "positive integer" but got "invalid format '%d'; use format %f, %e, %g or %a for numeric objects"
- **Root Cause**: Input validation error messages have changed - tests expect different error messages

#### Test: "apply_intrinsic_identifiability_core zeros tiny HRFs" (Line 779)
- **Error**: `result$Xi_ident_matrix[, tiny_voxel] not equal to rep(0, m)`
- **Differences**: 3/3 mismatches with average diff of 0.215
- **Root Cause**: Tiny HRF values (1e-12) are not being properly zeroed out as expected

## Summary of Issues

1. **Numerical Stability**: Several tests involve singular or near-singular matrices that cause failures in matrix decompositions (chol, qr.solve)

2. **Error Message Mismatches**: Some tests are checking for specific error messages that have changed in the implementation

3. **Numerical Precision**: The identifiability constraints are not meeting the expected precision tolerances in some cases

4. **Edge Case Handling**: Tests specifically designed to test edge cases (rank deficient matrices, tiny values) are revealing that the implementation doesn't handle these cases as expected

## Recommendations

1. Add more robust handling of singular matrices in the LSS and gamma solving routines
2. Update test expectations to match current error messages
3. Review numerical tolerances in identifiability constraints
4. Consider adding regularization or fallback strategies for edge cases