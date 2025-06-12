# Testing Priorities for manifoldhrf Package

Based on architectural analysis with Gemini, this document outlines the testing priorities for the manifoldhrf package.

## Algorithmic Importance Ranking

Files ranked by their criticality to the method's correctness:

1. **`core_manifold_construction.R`** - The heart of the package, defines the HRF manifold using diffusion maps
2. **`core_voxelfit_engine.R`** - Projects voxel data onto manifold basis using SVD/GLM
3. **`core_spatial_smoothing.R`** - Graph Laplacian regularization (not just smoothing, but part of the model)
4. **`core_alternating_optimization.R`** - Orchestrates iterative refinement between components
5. **`core_lss.R`** - Generates trial-wise estimates that feed the manifold engine

## Testing Priority Ranking

Files ranked by where bugs are most likely or would have catastrophic consequences:

1. **`input_validation.R`** - The gatekeeper; prevents cascading errors
2. **`core_manifold_construction.R`** - Complex math requiring numerical precision
3. **`mhrf_lss.R`** - Integration point where data flow errors occur
4. **`core_voxelfit_engine.R`** - Heavy numerical lifting prone to instability
5. **`core_lss.R`** - Complex trial iteration logic prone to off-by-one errors

## Critical Test Scenarios

### 1. Input Validation (`input_validation.R`)

**Test: Mismatched Dimensions**
- Scenario: 4D BOLD with mask of different size, event file with wrong number of timepoints
- Failure: Misaligned data producing nonsense results that look plausible
- Real trigger: Manual reslicing, aborted acquisitions

**Test: Pathological Values**
- Scenario: NA, NaN, Inf in voxel time series
- Failure: Silent propagation corrupting entire analysis
- Real trigger: De-spiking artifacts, sampling outside brain

**Test: Inconsistent Factor Levels**
- Scenario: 'face' vs 'Face', trailing whitespace, missing conditions
- Failure: Split conditions reducing power, rank-deficient designs
- Real trigger: Manual CSV editing, spreadsheet auto-formatting

### 2. Manifold Construction (`core_manifold_construction.R`)

**Test: Degenerate Data**
- Scenario: All voxels perfectly correlated, constant trial responses
- Failure: Zero distances, eigendecomposition crashes, collapsed manifold
- Real trigger: Failed LSS leaving near-zero betas, single tissue sampling

**Test: Disconnected Components**
- Scenario: Two well-separated clusters in beta space
- Failure: Inf geodesic distances, artificial bridges distorting structure
- Real trigger: Alert vs drowsy states, high within-block similarity

**Test: High Dimensional Curse**
- Scenario: 40,000 voxels but only 200 trials
- Failure: Unstable manifold, pure noise interpreted as structure
- Real trigger: Default neuroimaging state - must have robust dimensionality reduction

### 3. Pipeline Integration (`mhrf_lss.R`)

**Test: LSS Failure Propagation**
- Scenario: Rank-deficient trials producing NA/zero betas
- Failure: Passing bad data to manifold causing crashes or distortion
- Real trigger: Motion spikes perfectly explaining trial variance

**Test: Inconsistent Masking**
- Scenario: Different masks between LSS and manifold stages
- Failure: Dimension mismatches, scrambled voxel identities
- Real trigger: Whole-brain GLM but ROI-specific manifold

**Test: Metadata Alignment**
- Scenario: Dropped trials creating non-sequential numbering
- Failure: Misaligned brain-behavior relationships
- Real trigger: Behavioral data cleaning removing invalid trials

### 4. Voxel Fit Engine (`core_voxelfit_engine.R`)

**Test: Perfect Collinearity**
- Scenario: Design column is exact linear combination of others
- Failure: "System is singular" crash user can't debug
- Real trigger: Redundant regressors, dependent tissue signals

**Test: Near Collinearity**
- Scenario: Correlations of 0.999 between regressors
- Failure: Wildly unstable betas (±1e9), huge standard errors
- Real trigger: Overlapping event-related responses

**Test: Zero Variance Voxels**
- Scenario: Constant time series
- Failure: Division by zero in t-stats, Inf/NaN in maps
- Real trigger: Voxels outside brain, imperfect masks

### 5. LSS Implementation (`core_lss.R`)

**Test: Rapid Overlapping Events**
- Scenario: Trials every 2s with 16s HRF
- Failure: High collinearity systematically underestimating amplitudes
- Real trigger: Any rapid event-related design

**Test: Single-Trial Conditions**
- Scenario: Oddball with one instance
- Failure: Undefined nuisance regressor, misspecified model
- Real trigger: Rare events, error trials, target detection

**Test: Boundary Events**
- Scenario: Events at first/last TRs
- Failure: Truncated HRF causing systematic underestimation
- Real trigger: Designs without sufficient padding

## Consolidation Opportunities

1. **Voxel fitting files** (`core_voxelfit_engine.R`, `core_voxelwise_fit.R`, `robust_voxelwise_fit.R`) - Could be unified with method parameter
2. **`soundness_improvements.R`** - Should be eliminated, contents distributed to appropriate modules
3. **Interface files** - Could potentially merge `mhrf_lss_interface.R` and `simple_interface.R`

## Implementation Priority

1. Create comprehensive test suite for `input_validation.R` covering all edge cases
2. Build numerical precision tests for `core_manifold_construction.R`
3. Develop integration tests for `mhrf_lss.R` data flow
4. Add stability tests for `core_voxelfit_engine.R`
5. Create edge case tests for `core_lss.R` event handling