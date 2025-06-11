# LSS Memory Optimization Summary

## 1. Concrete Implementation Example

The key modification to `run_lss_voxel_loop_core` is to actually use the `ram_heuristic_GB_for_Rt` parameter to decide between three strategies:

```r
# Calculate memory requirement
precompute_memory_gb <- (n * T_trials * V * 8) / (1024^3)

# Decision logic
if (precompute_memory_gb <= ram_heuristic_GB_for_Rt) {
  # STRATEGY 1: Full precomputation (current implementation)
  all_C <- array(0, dim = c(n, T_trials, V))
  # ... existing code
} else if (optimal_chunk_size >= 100) {
  # STRATEGY 2: Chunked processing
  # Process voxels in memory-efficient chunks
} else {
  # STRATEGY 3: Pure streaming
  # Compute convolved regressors on-the-fly for each voxel
}
```

See `R/core_lss_improved.R` for the full implementation.

## 2. Heuristic for Deciding When to Precompute vs Stream

The recommended heuristic considers:

### A. Dataset Size vs Available Memory
```r
data_gb <- (n_timepoints * n_trials * n_voxels * 8) / (1024^3)
system_gb <- get_system_memory()

if (data_gb < system_gb * 0.25) {
  # Precompute - fits comfortably
} else if (data_gb < system_gb * 0.5) {
  # Chunk - moderate fit
} else {
  # Stream - too large
}
```

### B. Practical Thresholds
- **Small datasets (< 1 GB)**: Always precompute
- **Medium datasets (1-10 GB)**: Use chunking with 500-1000 voxels per chunk
- **Large datasets (> 10 GB)**: Stream computation

### C. Performance vs Memory Trade-off
- **Precomputation**: ~3-5x faster but requires T×V×n×8 bytes
- **Chunking**: ~1.5-2x slower but limits memory to chunk_size×T×n×8 bytes
- **Streaming**: ~5-10x slower but minimal memory overhead

## 3. Memory-Efficient Alternatives to 3D Array

Instead of creating `all_C <- array(0, dim = c(n, T_trials, V))`:

### A. Chunked Processing (Recommended)
```r
# Process voxels in chunks of 500-1000
chunk_size <- min(1000, floor(ram_limit_gb * 1e9 / (n * T_trials * 8)))
for (chunk in chunks) {
  chunk_C <- array(0, dim = c(n, T_trials, length(chunk)))
  # Process chunk
  rm(chunk_C); gc()
}
```

### B. Sparse Matrix Storage
```r
# When trial matrices are sparse (common in event-related designs)
sparse_C_list <- lapply(1:T_trials, function(t) {
  X_sparse <- Matrix::Matrix(X_trials[[t]], sparse = TRUE)
  C_t <- X_sparse %*% H_matrix
  return(Matrix::Matrix(C_t, sparse = TRUE))
})
```

### C. Memory-Mapped Files (for very large datasets)
```r
# Use disk-backed arrays
if (requireNamespace("ff", quietly = TRUE)) {
  all_C <- ff::ff(0, dim = c(n, T_trials, V), vmode = "double")
  # Process normally, but data lives on disk
}
```

## 4. Using Future Package for Memory-Aware Scheduling

The future package can automatically manage memory limits:

```r
# Set up memory-limited workers
future::plan(
  future::multisession,
  workers = n_jobs,
  globals = list(
    maxSize = 2 * 1024^3  # 2 GB per worker
  )
)

# Future will automatically schedule work to respect memory limits
voxel_futures <- lapply(1:V, function(v) {
  future::future({
    # Voxel processing code
  }, globals = list(size = estimated_size))
})
```

## Implementation Recommendations

1. **Update `run_lss_voxel_loop_core` to use the decision logic** shown in `core_lss_improved.R`

2. **Add an "auto" option** for `ram_heuristic_GB_for_Rt`:
   ```r
   if (ram_heuristic_GB_for_Rt == "auto") {
     ram_heuristic_GB_for_Rt <- calculate_optimal_ram_heuristic(...)
   }
   ```

3. **Provide presets** for common scenarios:
   - `"fast"`: Maximize speed (use all available RAM)
   - `"balanced"`: Good performance with safety margin
   - `"memory"`: Minimize memory usage

4. **Add progress reporting** that shows which strategy was chosen

5. **Consider adding a `max_chunk_size` parameter** to give users more control

The key insight is that the current implementation always does full precomputation, but the `ram_heuristic_GB_for_Rt` parameter exists specifically to enable smarter memory management - it just needs to be implemented!