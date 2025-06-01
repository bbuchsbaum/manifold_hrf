# Diagnostic Tests Summary

## Overview
Added two comprehensive diagnostic tests for the M-HRF-LSS core algorithm that validate mathematical correctness and algorithmic stability.

## Test 1: Signal Reconstruction Fidelity and Manifold Geometry

This test verifies:
1. **Manifold Embedding Quality**: The diffusion map embedding preserves HRF relationships
2. **Reconstruction Accuracy**: HRFs can be accurately reconstructed from manifold coordinates
3. **Noise Robustness**: The method is stable under different SNR conditions

### Key Validations:
- Eigenvalue decay of the manifold (smooth decay indicates good manifold)
- HRF reconstruction error (< 10% Frobenius norm error)
- Neighborhood preservation (> 70% of k-nearest neighbors preserved)
- Recovery quality at multiple SNR levels (Inf, 2, 1, 0.5)
- Monotonic degradation with decreasing SNR

### Test Design:
- Creates a structured HRF library with systematic variations in peak time and undershoot
- Tests manifold construction and reconstruction pipeline
- Evaluates performance under controlled noise conditions

## Test 2: Trial-wise Estimation Accuracy and Efficiency

This test verifies:
1. **Unbiased Estimation**: Trial-wise beta estimates are unbiased
2. **Numerical Equivalence**: Woodbury implementation matches naive LSS
3. **Efficiency**: Computational gains without sacrificing accuracy
4. **Overlap Handling**: Correctly handles temporally overlapping trials

### Key Validations:
- Woodbury vs naive LSS agreement (< 1e-10 difference)
- Beta estimate bias (< 0.05 mean error)
- Correlation with ground truth (> 0.85 at SNR=1.5)
- Variance explained (> 70%)
- Condition-level aggregation accuracy (> 0.95 correlation)
- Variance inflation bounds (< 3x for realistic scenarios)
- Overlapping trial discrimination

### Test Design:
- Simulates realistic fMRI experiment with overlapping trials
- Compares Woodbury matrix identity implementation against naive approach
- Tests recovery of known trial-wise amplitudes
- Evaluates handling of confounds and noise

## Current Status

The tests are partially passing but reveal some areas for improvement:

### Issues Identified:
1. **Manifold reconstruction error**: Higher than expected (89% vs target 90%)
2. **Index bounds error**: Issue with manifold coordinate distance computation
3. **Woodbury precision**: Numerical differences larger than expected
4. **Recovery performance**: Lower correlation than target at SNR=1.5
5. **Variance inflation**: Some voxels exceed the 3x threshold

### Strengths Confirmed:
- Basic algorithm structure is sound
- Manifold construction works
- LSS implementation functions
- Handles multiple conditions and trials

## Recommendations

1. **Improve numerical stability** in the Woodbury implementation
2. **Add regularization** to the manifold construction for better conditioning
3. **Tune parameters** for the specific test scenarios
4. **Consider adaptive thresholds** based on data characteristics

## Value of These Tests

These diagnostic tests provide:
1. **Ground truth validation** with known signals
2. **Mathematical correctness** verification
3. **Performance benchmarks** for optimization
4. **Edge case identification** for robustness improvements
5. **Quantitative metrics** for algorithm comparison

The tests serve as both validation tools and benchmarks for future improvements to the M-HRF-LSS algorithm.