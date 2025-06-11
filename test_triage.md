# Test Failure Triage

## Category 1: Ghosts of Woodbury (DELETE THESE)

1. **test-core-algorithm-diagnostics.R:416** - Calls `run_lss_woodbury_corrected` directly
2. **test-debug-woodbury.R:55** - Entire file is about Woodbury debugging
3. **test-lss-correction-validation.R:206** - Calls `run_lss_woodbury_corrected` directly
4. **test-fmrireg-benchmarks.R:397** - Calls `run_lss_for_voxel_corrected` (removed alias)

## Category 2: API Mismatches (UPDATE THESE)

5. **test-core-voxelfit-engine.R:71** - Calls `prepare_projection_matrix` (need to update)
6. **test-lss-correction-validation.R:50** - Calls `prepare_projection_matrix` (need to update)
7. **test-lss-formula-verification.R:88** - Calls `run_lss_for_voxel_corrected_full` (removed)
8. **test-lss-loop-core.R:29** - Calls `prepare_projection_matrix`
9. **test-lss-loop-core.R:87** - Calls `prepare_projection_matrix`
10. **test-lss-lsa-equivalence.R:68** - Wrong argument `Z_confounds` to `run_lss_for_voxel`

## Category 3: Legitimate Failures (DEBUG THESE)
None in the first 10 - they're all missing function errors!