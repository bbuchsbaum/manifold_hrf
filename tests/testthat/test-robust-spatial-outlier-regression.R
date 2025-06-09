library(testthat)

# helper functions replicating previous voxel-wise implementations

detect_outlier_timepoints_loop <- function(Y_data, threshold = 3, min_weight = 0.1) {
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  weights <- matrix(1, n, V)
  for (v in seq_len(V)) {
    y <- Y_data[, v]
    y_med <- median(y)
    y_mad <- median(abs(y - y_med)) * 1.4826
    if (y_mad > .Machine$double.eps) {
      z <- abs(y - y_med) / y_mad
      idx <- which(z > threshold)
      if (length(idx) > 0) {
        weights[idx, v] <- pmax(min_weight, 1 - (z[idx] - threshold) / threshold)
      }
    }
  }
  weights
}

screen_voxels_loop <- function(Y_data, min_variance = 1e-6, max_spike_fraction = 0.1) {
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  keep_voxel <- rep(TRUE, V)
  flag_voxel <- rep(FALSE, V)
  quality_scores <- numeric(V)
  for (v in seq_len(V)) {
    y <- Y_data[, v]
    if (any(!is.finite(y))) {
      keep_voxel[v] <- FALSE
      quality_scores[v] <- 0
      flag_voxel[v] <- TRUE
      next
    }
    y_var <- var(y)
    if (is.na(y_var) || y_var < min_variance) {
      keep_voxel[v] <- FALSE
      quality_scores[v] <- 0
      next
    }
    y_diff <- abs(diff(y))
    y_diff_med <- median(y_diff)
    y_diff_mad <- median(abs(y_diff - y_diff_med)) * 1.4826
    if (y_diff_mad > 0) {
      spike_threshold <- y_diff_med + 5 * y_diff_mad
      n_spikes <- sum(y_diff > spike_threshold)
      spike_fraction <- n_spikes / (n - 1)
      if (spike_fraction > max_spike_fraction) {
        flag_voxel[v] <- TRUE
      }
    }
    quality_scores[v] <- 1 / (1 + mean(y_diff^2) / y_var)
  }
  list(
    keep = keep_voxel,
    flag = flag_voxel,
    quality_scores = quality_scores,
    n_excluded = sum(!keep_voxel),
    n_flagged = sum(flag_voxel)
  )
}


test_that("detect_outlier_timepoints matches previous implementation", {
  set.seed(123)
  Y <- matrix(rnorm(40), 10, 4)
  Y[3, 1] <- Y[3, 1] + 5
  Y[7, 3] <- Y[7, 3] - 5
  ref <- detect_outlier_timepoints_loop(Y, threshold = 3)
  res <- suppressMessages(detect_outlier_timepoints(Y, threshold = 3))
  expect_equal(res, ref)
})


test_that("screen_voxels matches previous implementation", {
  set.seed(42)
  Y <- matrix(rnorm(60), 20, 3)
  Y[,1] <- 0
  Y[5,2] <- Inf
  Y[15,3] <- -Inf
  Y[10,3] <- Y[10,3] + 10
  ref <- suppressWarnings(screen_voxels_loop(Y))
  res <- suppressWarnings(screen_voxels(Y))
  expect_equal(res$keep, ref$keep)
  expect_equal(res$flag, ref$flag)
  expect_equal(res$quality_scores, ref$quality_scores)
  expect_equal(res$n_excluded, ref$n_excluded)
  expect_equal(res$n_flagged, ref$n_flagged)
})

