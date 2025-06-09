context("compute_local_snr")

test_that("vectorized compute_local_snr matches loop version", {
  set.seed(123)
  n <- 40
  V <- 8
  Y <- matrix(rnorm(n * V), n, V)
  Ypred <- matrix(rnorm(n * V), n, V)

  loop_temporal <- function(Y_data) {
    V <- ncol(Y_data)
    out <- numeric(V)
    for (v in seq_len(V)) {
      y <- Y_data[, v]
      signal_var <- var(y)
      noise_mad <- median(abs(diff(y))) * 1.4826
      noise_var <- noise_mad^2
      out[v] <- signal_var / (noise_var + .Machine$double.eps)
    }
    out[is.na(out)] <- 1
    out[is.infinite(out)] <- 100
    out[out > 100] <- 100
    out[out < 0.1] <- 0.1
    out
  }

  loop_residual <- function(Y_data, Y_predicted) {
    V <- ncol(Y_data)
    out <- numeric(V)
    for (v in seq_len(V)) {
      signal_var <- var(Y_predicted[, v])
      residuals <- Y_data[, v] - Y_predicted[, v]
      noise_var <- var(residuals)
      out[v] <- signal_var / (noise_var + .Machine$double.eps)
    }
    out[is.na(out)] <- 1
    out[is.infinite(out)] <- 100
    out[out > 100] <- 100
    out[out < 0.1] <- 0.1
    out
  }

  expect_equal(manifoldhrf:::compute_local_snr(Y),
               loop_temporal(Y))
  expect_equal(manifoldhrf:::compute_local_snr(Y, Ypred, method = "residual"),
               loop_residual(Y, Ypred))
})
