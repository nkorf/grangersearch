# Tests for granger_causality_test main functionality

test_that("granger_causality_test returns correct class", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result <- granger_causality_test(x = x, y = y)

  expect_s3_class(result, "granger_result")
})

test_that("granger_causality_test detects known causal relationship", {
  set.seed(42)
  n <- 200

  # Create X as random walk
  x <- cumsum(rnorm(n))

  # Y depends on lagged X (X causes Y)
  y <- c(0, 0.8 * x[1:(n-1)]) + rnorm(n, sd = 0.3)

  result <- granger_causality_test(x = x, y = y)

  # X should Granger-cause Y

  expect_true(result$x_causes_y)
  expect_lt(result$p_value_xy, 0.05)
})

test_that("granger_causality_test detects no causality for independent series", {
  set.seed(123)
  n <- 200

  # Two independent random walks
  x <- cumsum(rnorm(n))
  y <- cumsum(rnorm(n))

  result <- granger_causality_test(x = x, y = y)

  # Neither should cause the other (most of the time)
  # Using alpha = 0.01 for stricter test
  result_strict <- granger_causality_test(x = x, y = y, alpha = 0.01)
  expect_false(result_strict$x_causes_y && result_strict$y_causes_x)
})

test_that("granger_causality_test returns all expected components", {
  set.seed(123)
  x <- cumsum(rnorm(50))
  y <- c(0, x[1:49]) + rnorm(50, sd = 0.5)

  result <- granger_causality_test(x = x, y = y)

  expect_true("x_causes_y" %in% names(result))
  expect_true("y_causes_x" %in% names(result))
  expect_true("p_value_xy" %in% names(result))
  expect_true("p_value_yx" %in% names(result))
  expect_true("test_statistic_xy" %in% names(result))
  expect_true("test_statistic_yx" %in% names(result))
  expect_true("lag" %in% names(result))
  expect_true("alpha" %in% names(result))
  expect_true("test" %in% names(result))
  expect_true("n" %in% names(result))
  expect_true("x_name" %in% names(result))
  expect_true("y_name" %in% names(result))
  expect_true("call" %in% names(result))
})

test_that("granger_causality_test respects lag parameter", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result1 <- granger_causality_test(x = x, y = y, lag = 1)
  result2 <- granger_causality_test(x = x, y = y, lag = 3)

  expect_equal(result1$lag, 1L)
  expect_equal(result2$lag, 3L)
})

test_that("granger_causality_test respects alpha parameter", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result1 <- granger_causality_test(x = x, y = y, alpha = 0.05)
  result2 <- granger_causality_test(x = x, y = y, alpha = 0.01)

  expect_equal(result1$alpha, 0.05)
  expect_equal(result2$alpha, 0.01)

  # Same p-value, but different conclusions possible
  expect_equal(result1$p_value_xy, result2$p_value_xy)
})

test_that("granger_causality_test works with tidyverse syntax", {
  skip_if_not_installed("tibble")

  set.seed(123)
  df <- tibble::tibble(
    price = cumsum(rnorm(100)),
    volume = c(0, cumsum(rnorm(99)))
  )

  result <- granger_causality_test(df, price, volume)

  expect_s3_class(result, "granger_result")
  expect_equal(result$x_name, "price")
  expect_equal(result$y_name, "volume")
})

test_that("granger_causality_test works with pipe operator", {
  skip_if_not_installed("tibble")

  set.seed(123)
  df <- tibble::tibble(
    a = cumsum(rnorm(100)),
    b = c(0, cumsum(rnorm(99)))
  )

  # Using base pipe
  result <- df |> granger_causality_test(a, b)

  expect_s3_class(result, "granger_result")
  expect_equal(result$x_name, "a")
  expect_equal(result$y_name, "b")
})

test_that("tidy method returns correct tibble structure", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result <- granger_causality_test(x = x, y = y)
  tidy_result <- tidy(result)

  expect_s3_class(tidy_result, "tbl_df")
  expect_equal(nrow(tidy_result), 2)
  expect_true(all(c("direction", "cause", "effect", "statistic", "p.value", "significant") %in% names(tidy_result)))
})

test_that("glance method returns correct tibble structure", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result <- granger_causality_test(x = x, y = y)
  glance_result <- glance(result)

  expect_s3_class(glance_result, "tbl_df")
  expect_equal(nrow(glance_result), 1)
  expect_true(all(c("nobs", "lag", "alpha", "test", "bidirectional", "x_causes_y", "y_causes_x") %in% names(glance_result)))
})

test_that("p-values are between 0 and 1", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y)

  expect_gte(result$p_value_xy, 0)
  expect_lte(result$p_value_xy, 1)
  expect_gte(result$p_value_yx, 0)
  expect_lte(result$p_value_yx, 1)
})

test_that("test statistics are numeric", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y)

  expect_type(result$test_statistic_xy, "double")
  expect_type(result$test_statistic_yx, "double")
  expect_false(is.na(result$test_statistic_xy))
  expect_false(is.na(result$test_statistic_yx))
})
