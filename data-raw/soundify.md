# Soundify M-HRF-LSS: Low-Cost Robustness Improvements

## Goal
Make the M-HRF-LSS method more robust and increase likelihood of "just working" on real data without major architectural changes.

## Guiding Principles
- Minimal code changes
- Focus on numerical stability
- Add safety nets and fallbacks
- Improve default parameter selection
- Make failures graceful and informative

---

## Critical Soundness Improvements

### 1. **SOUND-MANIFOLD-STABILITY: Manifold Construction Robustness**
**Problem**: Manifold construction can fail with degenerate HRF libraries or produce unstable embeddings.

**Solutions**:
- Add minimum diversity check for HRF library (condition number check)
- Implement fallback to PCA if diffusion map fails or is ill-conditioned
- Add automatic duplicate HRF removal with warning
- Check for isolated points in affinity matrix and handle gracefully
- Add manifold quality metric (e.g., residual variance) to outputs

**Implementation**: ~20 lines in `calculate_manifold_affinity_core()` and `get_manifold_basis_reconstructor_core()`

---

### 2. **SOUND-PARAM-ADAPTIVE: Adaptive Default Parameters**
**Problem**: Fixed defaults may not work across different data characteristics.

**Solutions**:
- Scale `lambda_gamma` based on data variance: `lambda_gamma = 0.01 * var(Y_data)`
- Set `k_local_nn_for_sigma` as `max(7, round(sqrt(N_library)))`
- Adapt `lambda_spatial_smooth` based on voxel density: fewer voxels → less smoothing
- Add data-driven parameter suggestions to QC report
- Create `suggest_parameters()` function that analyzes data characteristics

**Implementation**: New helper function (~50 lines) + small modifications to defaults

---

### 3. **SOUND-VOXFIT-REGULARIZE: Improve Voxel-wise Fitting Stability**
**Problem**: SVD-based extraction can be unstable for poorly conditioned Gamma matrices.

**Solutions**:
- Add condition number check before SVD with automatic regularization increase
- Implement "soft" SVD truncation based on singular value gap
- Add fallback to first PC when SVD fails
- Track and report per-voxel fitting stability metric
- Option to use robust SVD (e.g., randomized SVD for large problems)

**Implementation**: ~30 lines in `extract_xi_beta_raw_svd_core()`

---

### 4. **SOUND-SPATIAL-ADAPTIVE: Adaptive Spatial Smoothing**
**Problem**: Fixed spatial smoothing may over/under-smooth depending on brain region.

**Solutions**:
- Make smoothing strength proportional to local SNR (low SNR → more smoothing)
- Implement edge-preserving smoothing (reduce smoothing at tissue boundaries)
- Add minimum neighbor check - skip smoothing for isolated voxels
- Allow spatially varying lambda based on local data quality
- Add diagnostic showing smoothing strength map

**Implementation**: ~40 lines enhancement to `apply_spatial_smoothing_core()`

---

### 5. **SOUND-HRF-CONSTRAINTS: Physiological HRF Constraints**
**Problem**: Optimization might produce non-physiological HRF shapes.

**Solutions**:
- Add soft constraints to keep HRFs positive after peak
- Penalize HRFs with peaks too early (<2s) or late (>10s)
- Ensure HRF integral is positive (net positive response)
- Add HRF "reasonableness" score to QC metrics
- Option to project HRFs back to physiological subspace

**Implementation**: ~30 lines in `apply_intrinsic_identifiability_core()` + new helper

---

### 6. **SOUND-INIT-SMART: Smart Initialization**
**Problem**: Poor initialization can lead to local minima.

**Solutions**:
- Initialize with standard GLM+canonical HRF for warm start
- Use spatial clustering to identify similar voxels for initialization
- Add option to initialize from previous run (for iterative refinement)
- Detect and reinitialize stuck voxels (very low variance in Xi)
- Report initialization quality metrics

**Implementation**: New `smart_initialize()` function (~60 lines) + integration hooks

---

### 7. **SOUND-OUTLIER-ROBUST: Outlier Detection and Handling**
**Problem**: Outlier voxels/timepoints can corrupt estimates.

**Solutions**:
- Add automatic outlier timepoint detection (>3 SD) with soft downweighting
- Implement voxel screening: skip voxels with too low temporal variance
- Add robust regression option (Huber weights) for high-noise data
- Flag and report outlier voxels/timepoints in QC
- Option to use median instead of mean in key computations

**Implementation**: ~40 lines across multiple functions + new helper

---

### 8. **SOUND-CONVERGE-CHECK: Convergence Diagnostics**
**Problem**: May not converge or converge to poor solutions.

**Solutions**:
- Add convergence tracking for Component 1 (not just Component 4)
- Implement early stopping based on relative change
- Add "solution quality" metric (e.g., reconstruction error + smoothness)
- Warn if convergence is suspiciously fast (might be stuck)
- Option for multiple random restarts with best solution selection

**Implementation**: ~30 lines of tracking code + convergence helper

---

### 9. **SOUND-FALLBACK-CASCADE: Graceful Degradation**
**Problem**: Complete failure leaves users without results.

**Solutions**:
- Implement fallback cascade: Manifold → PCA → Canonical HRF
- Add "method used" indicator for each voxel
- If spatial smoothing fails, continue without it
- If LSS fails for some trials, fall back to condition-level
- Always return something interpretable with quality indicators

**Implementation**: ~50 lines of try-catch logic with fallback paths

---

### 10. **SOUND-MEMORY-SAFE: Memory Safety Improvements**
**Problem**: Large datasets may cause memory failures.

**Solutions**:
- Add memory requirement estimation before processing
- Implement automatic chunking for voxel loops
- Add option for disk-backed matrices (via ff/bigmemory packages)
- Clear large intermediate objects proactively
- Add memory usage tracking to diagnostics

**Implementation**: Memory helpers (~40 lines) + gc() calls

---

## Quick Wins (< 10 lines each)

### 11. **SOUND-SCALE-ROBUST: Input Scaling Robustness**
- Auto-scale Y_data if values are too large/small (>1e6 or <1e-6)
- Warn if data appears to be integer-valued (not percent signal change)

### 12. **SOUND-RANK-CHECK: Design Matrix Rank Checking**
- Check and warn if design matrices are rank-deficient
- Automatic column removal for perfect collinearity

### 13. **SOUND-ZERO-HANDLE: Better Zero Handling**
- Check for all-zero voxels upfront and skip processing
- Handle zero-variance voxels gracefully

### 14. **SOUND-CONDITION-MONITOR: Condition Number Monitoring**
- Track condition numbers throughout pipeline
- Warn when approaching numerical limits

### 15. **SOUND-PROGRESS-FEEDBACK: Better Progress Feedback**
- Add progress bars for long operations
- Report estimated time remaining
- Allow graceful interruption with partial results

---

## Testing Improvements

### 16. **SOUND-TEST-ADVERSARIAL: Adversarial Test Cases**
- Test with pathological inputs (all zeros, all same value, single spike)
- Test with minimal data (3 voxels, 10 timepoints)
- Test with high collinearity designs
- Test with extreme noise levels

### 17. **SOUND-TEST-RECOVERY: Recovery Testing**
- Test that fallbacks actually work
- Test partial failure scenarios
- Test memory limits
- Test with missing data (NAs)

---

## Default Workflow Improvements

### 18. **SOUND-WORKFLOW-PRESET: Preset Configurations**
```r
# Add preset configurations
get_preset_params <- function(preset = c("conservative", "balanced", "aggressive")) {
  switch(preset,
    conservative = list(
      lambda_gamma = 0.1,
      lambda_spatial_smooth = 1.0,
      m_manifold_dim_target = 3
    ),
    balanced = list(
      lambda_gamma = 0.01,
      lambda_spatial_smooth = 0.5, 
      m_manifold_dim_target = 5
    ),
    aggressive = list(
      lambda_gamma = 0.001,
      lambda_spatial_smooth = 0.1,
      m_manifold_dim_target = 8
    )
  )
}
```

---

## Priority Order

**Immediate (High Impact, Low Effort):**
1. SOUND-MANIFOLD-STABILITY (prevents catastrophic failures)
2. SOUND-PARAM-ADAPTIVE (better out-of-box experience)
3. SOUND-VOXFIT-REGULARIZE (core stability)
4. SOUND-FALLBACK-CASCADE (always get results)

**Next Wave:**
5. SOUND-OUTLIER-ROBUST
6. SOUND-HRF-CONSTRAINTS
7. SOUND-SPATIAL-ADAPTIVE
8. SOUND-SCALE-ROBUST through SOUND-CONDITION-MONITOR

**Testing & Documentation:**
9. SOUND-TEST-ADVERSARIAL
10. SOUND-WORKFLOW-PRESET

---

## Success Metrics

After implementing these improvements, the method should:
- Never crash with real fMRI data
- Always produce interpretable results (even if degraded)
- Provide clear feedback about what happened
- Work reasonably well with default parameters
- Gracefully handle edge cases
- Give users confidence through transparency

## Estimated Total Effort
- ~400-500 lines of defensive code
- ~100 lines of new tests  
- 2-3 days of implementation
- Minimal architectural changes
- High impact on robustness