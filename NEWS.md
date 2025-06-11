# manifoldhrf 0.1.0.9000

## Breaking Changes

* **Complete LSS implementation overhaul**: All LSS (Least Squares Separate) functionality now routes through the external `fmrilss` package for improved performance and reliability.
* **Removed legacy functions**:
  - `run_lss_woodbury_corrected()` - Use `run_lss_for_voxel()` instead
  - `run_lsa_woodbury()` - Use `run_lss_for_voxel()` instead  
  - `prepare_projection_matrix()` - Pass confound matrix directly to LSS functions
  - All `run_lss_for_voxel_corrected*` aliases
  - Internal helper functions `.compute_lss_beta_series()`, `.lsa_beta()`, `.project_confounds()`
* **Simplified API**: The package now provides a cleaner, simpler interface for LSS computations

## Improvements

* Added version constraint for `fmrilss` dependency (>= 0.1.0)
* Updated documentation to reflect new fmrilss backend
* Improved numerical stability in tests by adding appropriate tolerances
* Removed obsolete test files (test-woodbury-mathematical-verification.R)

## Bug Fixes

* Fixed double projection issue when data is already projected
* Fixed Windows compatibility issue in parallel processing
* Fixed NaN values occurring when trials don't fit within time series

# manifoldhrf 0.1.0

* Initial release