# Tests for smart initialization and stuck voxel handling

test_that("smart_initialize provides reasonable starting values", {
  set.seed(123)
  n <- 60
  p <- 10
  V <- 6
  k <- 2
  m <- 3

  # canonical HRF
  t_hrf <- seq(0, (p - 1) * 0.5, by = 0.5)
  hrf <- dgamma(t_hrf, shape = 6, rate = 1)
  hrf <- hrf / sum(hrf)

  # simple event designs
  X_list <- lapply(seq_len(k), function(c) {
    X <- matrix(0, n, p)
    onsets <- seq(5 + c, by = 20, length.out = 3)
    for (o in onsets) {
      if (o + p - 1 <= n) {
        X[o:(o + p - 1), ] <- X[o:(o + p - 1), ] + diag(p)
      }
    }
    X
  })

  # generate data with varying SNR
  Y <- matrix(0, n, V)
  betas <- matrix(runif(k * V, 0.5, 1.5), k, V)
  noise_sd <- c(rep(0.1, V / 2), rep(1, V / 2))
  for (v in seq_len(V)) {
    for (c in seq_len(k)) {
      Y[, v] <- Y[, v] + X_list[[c]] %*% (hrf * betas[c, v])
    }
    Y[, v] <- Y[, v] + rnorm(n, sd = noise_sd[v])
  }

  init <- smart_initialize(
    Y_data = Y,
    X_condition_list = X_list,
    hrf_canonical = hrf,
    use_spatial_clusters = FALSE,
    m_manifold_dim = m
  )

  expect_equal(dim(init$Xi_init), c(m, V))
  expect_equal(dim(init$Beta_init), c(k, V))
  expect_length(init$R2_init, V)
  expect_true(all(init$R2_init >= 0 & init$R2_init <= 1))
  expect_true(all(init$good_voxels %in% seq_len(V)))

  mean_good <- mean(init$R2_init[init$good_voxels])
  mean_bad <- mean(init$R2_init[-init$good_voxels])
  expect_gt(mean_good, mean_bad)
})


test_that("fix_stuck_voxels reinitializes low-variance columns", {
  set.seed(456)
  m <- 3
  V <- 4
  Xi_prev <- matrix(rnorm(m * V), m, V)
  Xi_curr <- Xi_prev
  Xi_curr[, c(1, 3)] <- Xi_curr[, c(1, 3)] + matrix(rnorm(m * 2, sd = 0.1), m, 2)

  Xi_fixed <- fix_stuck_voxels(
    Xi_current = Xi_curr,
    Xi_previous = Xi_prev,
    variance_threshold = 1e-4,
    reinit_sd = 0.5
  )

  diff_stuck <- colSums((Xi_fixed[, c(2, 4)] - Xi_prev[, c(2, 4)])^2)
  diff_ok <- colSums((Xi_fixed[, c(1, 3)] - Xi_curr[, c(1, 3)])^2)

  expect_true(all(diff_stuck > 1e-4))
  expect_true(all(diff_ok < 1e-6))
  expect_equal(dim(Xi_fixed), dim(Xi_curr))
})
