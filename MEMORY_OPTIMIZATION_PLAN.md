# Memory Optimization Plan for manifoldhrf

Based on Gemini's analysis of memory limitations in the package.

## Key Finding: Major Memory Bottleneck

The LSS trial precomputation creates an array of size `n × T × V`:
- **Example**: 300 timepoints × 100 trials × 50,000 voxels × 8 bytes = **11.2 GB**
- This is often larger than available RAM on typical systems

## Current Problem

The `ram_heuristic_GB_for_Rt` parameter exists but isn't being used:
```r
# Current code always does full precomputation:
all_C <- array(0, dim = c(n, T_trials, V))
```

## Gemini's Solution: Three-Strategy Approach

### 1. Full Precomputation (When Memory Allows)
- Use when: `required_memory < ram_heuristic_GB_for_Rt`
- Performance: Fastest
- Memory: Highest (11.2 GB for example dataset)

### 2. Chunked Processing (Moderate Memory)
- Use when: `required_memory < ram_heuristic_GB_for_Rt × 10`
- Process 500-1000 voxels at a time
- Performance: Good
- Memory: ~1-2 GB

### 3. Streaming (Minimal Memory)
- Use when: Memory is severely limited
- Compute each voxel on-demand
- Performance: Slower but still parallelizable
- Memory: Minimal (~100 MB)

## Implementation Created

File: `R/core_lss_memory_optimized.R` contains:

1. **`calculate_optimal_ram_heuristic()`** - Intelligently sets memory limits based on:
   - Available system memory
   - Dataset size
   - User strategy (fast/balanced/memory)

2. **`run_lss_voxel_loop_memory_optimized()`** - Drop-in replacement that:
   - Actually uses the `ram_heuristic_GB_for_Rt` parameter
   - Automatically selects the best strategy
   - Provides progress feedback

3. **Future package integration** - Leverages existing dependency for memory-aware scheduling

## Other Memory Optimizations Recommended

### High Priority
1. **Sparse matrices for trial designs** - Many fMRI designs are sparse
2. **In-place operations** - Avoid unnecessary copies
3. **Memory warnings** - Alert users before large allocations

### Medium Priority
1. **HDF5 support** - For datasets larger than RAM
2. **Single precision option** - Halve memory usage where appropriate
3. **Smarter parallelization** - Adjust workers based on available memory

### Example Usage

```r
# Auto-select strategy based on available memory
result <- run_lss_voxel_loop_memory_optimized(
  Y_proj_matrix = Y_data,
  X_trial_onset_list_of_matrices = X_list,
  H_shapes_allvox_matrix = H_shapes,
  # ... other parameters ...
  ram_heuristic_GB_for_Rt = NULL  # Auto-calculate
)

# Or specify memory limit
result <- run_lss_voxel_loop_memory_optimized(
  # ... parameters ...
  ram_heuristic_GB_for_Rt = 4.0  # Use max 4 GB
)

# Or use future-based version
result <- run_lss_with_future_memory_limits(
  Y_proj_matrix = Y_data,
  X_trial_onset_list_of_matrices = X_list,
  H_shapes_allvox_matrix = H_shapes,
  max_memory_gb = 8.0
)
```

## Expected Impact

- **Memory usage**: 11.2 GB → 1-2 GB (chunked) or ~100 MB (streaming)
- **Performance**: Minimal impact for chunked approach
- **Scalability**: Can now handle datasets 10x larger
- **User experience**: Automatic strategy selection with informative messages

## Next Steps

1. Test the memory-optimized implementation with real datasets
2. Benchmark performance vs memory trade-offs
3. Integrate into main `mhrf_analyze()` pipeline
4. Add memory profiling to QC reports
5. Document memory management in vignettes

---
*Analysis and recommendations by Gemini*