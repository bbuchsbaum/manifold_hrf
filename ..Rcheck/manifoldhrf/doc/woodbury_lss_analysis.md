# Woodbury LSS Implementation Analysis

## Executive Summary

The current LSS implementation in `run_lss_for_voxel_core()` is implementing a specific variant from the M-HRF-LSS proposal, not the standard Woodbury matrix identity. While the implementation follows the proposal's formula, it produces results that differ significantly from standard LSS.

## Key Findings

### 1. Mathematical Formulation

The implementation uses:
```
alpha_v_row = (1 - pc_v_row) / cv_v_row
S_effective_regressors_v = V_regressors_v * alpha_v_row + p_lss_vec
beta = S_effective_regressors_v' * Y_proj
```

This is NOT the standard Woodbury matrix identity:
```
(A + uv')^(-1) = A^(-1) - A^(-1)uv'A^(-1) / (1 + v'A^(-1)u)
```

### 2. Test Results

Comparing different methods on the same data:
- **True betas**: [1, -0.5, 2, 0, 1.5, -1, 0.8, 0.3]
- **Direct LSS**: [0.903, -0.469, 1.9, 0.875, 1.462, -1.71, 0.663, 0]
- **Corrected**: [0.903, -0.469, 1.9, 0.875, 1.462, -1.71, 0.663, 0]
- **Current Woodbury**: [0.991, -0.958, 1.726, 0.56, 1.448, -1.997, 0.869, 0]

Differences from direct LSS:
- Corrected: ~1e-6 (essentially perfect)
- Current: 0.09-0.49 (significant errors)

### 3. Conceptual Issues

1. **LSS vs Simultaneous Estimation**: 
   - True LSS estimates each trial separately with others as nuisance
   - The corrected implementation does simultaneous estimation
   - The current implementation attempts something in between

2. **Numerical Stability**:
   - When `cv_v_row` is near zero, `alpha_v_row` becomes huge (4.5e+15 observed)
   - This causes numerical instability and incorrect results

3. **Data Flow Confusion**:
   - The implementation expects projected Y but unprojected X
   - This hybrid approach may be causing issues

## Recommendations

### Option 1: Fix Current Implementation
The formula might be correct but has implementation issues:
- Add numerical safeguards for near-zero `cv_v_row`
- Verify the mathematical derivation
- Add comprehensive unit tests

### Option 2: Replace with Standard Approach
Implement standard LSS using established linear algebra:
```r
# For each trial t:
X_full = [x_t, X_others, Z]
beta_t = (X'X + λI)^(-1) X'y
```

### Option 3: Hybrid Approach
Keep the proposal's approach but add fallback:
- Use the alpha formula when numerically stable
- Fall back to standard approach when unstable

## Mathematical Verification

The proposal's formula appears to be a computational shortcut, but it's not equivalent to standard LSS. The derivation needs to be checked:

1. Is it trying to avoid matrix inversions?
2. What assumptions does it make?
3. Under what conditions is it equivalent to standard LSS?

## Conclusion

The current implementation has mathematical issues that cause significant errors compared to standard LSS. While it follows the proposal's formula, that formula itself may need revision or the implementation needs numerical safeguards.

The 0.09-0.49 differences from ground truth are too large for production use. Either the mathematical approach needs to be fixed or replaced with a standard, well-tested method.