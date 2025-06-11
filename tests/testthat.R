# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(manifoldhrf)

# Configure backtrace display: keep for errors, suppress for warnings
# - Errors show simplified backtraces for debugging
# - Warnings suppress backtraces to reduce clutter
options(
  rlang_backtrace_on_error = "branch",          # Show simplified backtrace for errors
  rlang_backtrace_on_warning = "none",          # Suppress backtrace for warnings
  rlang_backtrace_on_error_report = "branch",   # Show simplified backtrace for errors in reports
  rlang_backtrace_on_warning_report = "none"    # Suppress backtrace for warnings in reports
)

test_check("manifoldhrf")