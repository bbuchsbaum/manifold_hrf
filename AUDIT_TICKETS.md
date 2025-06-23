# Manifold HRF - Critical Issue Tickets

Based on the audit and test results, here are the critical issues that need immediate attention:

## CRITICAL-1: Memory Overflow Risk in Distance Calculations
**Severity**: Critical  
**Files Affected**: 
- `config_helpers.R:31-33`
- `distance_fallbacks.R:12`
- `core_manifold_construction.R:174-231`

**Issue**: Creating full distance matrices without memory checks can crash the system
**Current Behavior**: `as.matrix(dist(voxel_coords))` creates V×V matrix regardless of size
**Fix Required**:
```r
# Add memory check before distance calculation
check_distance_memory <- function(n_points) {
  estimated_gb <- (n_points * n_points * 8) / 1e9
  available_gb <- as.numeric(memory.limit()) / 1024
  if (estimated_gb > available_gb * 0.5) {
    stop(sprintf("Distance matrix would require %.1f GB (available: %.1f GB)", 
                 estimated_gb, available_gb))
  }
}
```

---

## CRITICAL-2: Numerical Instability in Manifold Construction
**Severity**: Critical  
**File**: `core_manifold_construction.R:142-160`
**Issue**: Division by near-zero sigma values causing NaN/Inf propagation
**Test Status**: Tests pass but with warnings about negative eigenvalues
**Fix Required**:
```r
# Line 142-160: Add minimum threshold
sigma_i[i] <- if (length(nz) > 0) max(median(nz), 1e-6) else 1e-6
```

---

## HIGH-1: Unsafe Variance Calculations
**Severity**: High  
**File**: `input_validation.R:79-80`
**Issue**: `var()` can return NaN for constant columns even with `na.rm=TRUE`
**Fix Required**:
```r
voxel_vars <- apply(Y_matrix, 2, function(x) {
  finite_vals <- x[is.finite(x)]
  if (length(finite_vals) < 2) return(0)
  var(finite_vals)
})
```

---

## HIGH-2: Missing Input Validation in C++ Functions
**Severity**: High  
**Files**: `distance_knn.cpp`, `RcppExports.R`
**Issue**: No validation for empty matrices or dimension mismatches
**Note**: C++ functions exist (contrary to audit finding) but lack input validation
**Fix Required**: Add checks in C++ code for:
- Empty matrices
- k > number of data points
- Dimension mismatches

---

## HIGH-3: Thread Safety in Logger
**Severity**: High  
**File**: `logger.R:11-26`
**Issue**: Logger not thread-safe but package uses parallel processing
**Fix Required**: Document thread safety limitations or implement mutex

---

## MEDIUM-1: API Naming Clarity
**Severity**: Medium  
**Files**: `mhrf_lss_interface.R:76`, `simple_interface.R:36`
**Issue**: Confusing function names (`mhrf_lss` vs `mhrf_lss_core`)
**Current State**: Actually not a conflict - functions have different names
**Recommendation**: Document the interface hierarchy clearly

---

## Test Results Summary

### Passing Tests:
- ✅ Core manifold construction (43 tests, 9 warnings)
- ✅ Input validation (9 tests, 0 warnings)

### Test Warnings:
- Negative eigenvalues in row-stochastic matrices (expected behavior)
- Manual sparsification inefficiency warnings
- LSS memory strategy tests generate 10,000+ warnings (needs investigation)

### Overall Status:
- Core functionality appears to work correctly
- Many warnings indicate suboptimal implementations
- Critical memory and numerical stability issues need fixes before production

## Recommended Action Plan

1. **Immediate** (1 day):
   - Add memory overflow protection (CRITICAL-1)
   - Fix numerical stability in sigma calculations (CRITICAL-2)
   - Fix variance calculation safety (HIGH-1)

2. **Short-term** (2-3 days):
   - Add C++ input validation (HIGH-2)
   - Document thread safety limitations (HIGH-3)
   - Investigate excessive warnings in LSS tests

3. **Medium-term** (1 week):
   - Refactor sparse matrix handling to avoid warnings
   - Improve documentation of interface hierarchy
   - Add integration tests for memory limits