# Updated Test Failure Triage

## Fixed Tests ✓
1. test-lss-correction-validation.R - Updated to use fmrilss directly
2. test-core-voxelfit-engine.R - Updated to use new interface

## Remaining Failures (17)

### Files to Update or Delete:
1. **test-lss-separate.R** - 2 failures
   - Line 51: run_lss_woodbury_corrected matches LS-S
   - Line 91: run_lss_for_voxel_corrected_full matches LS-S
   
2. **test-lss.R** - 6 failures
   - Line 231: run_lss_for_voxel_corrected_full works correctly
   - Line 265: run_lss_for_voxel_corrected_full validates inputs
   - Line 339: run_lss_voxel_loop_corrected_test produces consistent results
   - Line 364: run_lss_voxel_loop_corrected_test validates inputs
   - Line 468: LSS integration test with known signal
   - Line 515: rank deficient trial regressors trigger warning

3. **test-lss-formula-verification.R** - 1 failure
   - Line 88: Current LSS formula matches reference implementation

4. **test-lss-loop-core.R** - 2 failures
   - Line 29: run_lss_voxel_loop_core matches single voxel implementation
   - Line 87: run_lss_voxel_loop_core works without precomputation

5. **test-lss-lsa-equivalence.R** - 1 failure
   - Line 68: LS-A and LS-S are equivalent for T=2 with confounds

6. **test-fmrireg-benchmarks.R** - 1 failure (duplicate)
   - Line 417: M-HRF-LSS performs trial-wise estimation correctly

## Action Plan:
1. Delete test-lss-separate.R (tests removed implementation)
2. Update test-lss.R to use new interface
3. Update or delete test-lss-formula-verification.R
4. Update test-lss-loop-core.R to test current implementation
5. Update test-lss-lsa-equivalence.R 
6. Update test-fmrireg-benchmarks.R line 417