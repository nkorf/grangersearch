# Tests for input validation

test_that("error on non-numeric x", {
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = letters[1:100], y = y),
    "must be numeric"
  )

  expect_error(
    granger_causality_test(x = as.factor(1:100), y = y),
    "must be numeric"
  )
})

test_that("error on non-numeric y", {
  x <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = letters[1:100]),
    "must be numeric"
  )
})

test_that("error on different length vectors", {
  x <- rnorm(100)
  y <- rnorm(50)

  expect_error(
    granger_causality_test(x = x, y = y),
    "must have the same length"
  )
})

test_that("error on too-short vectors", {
  x <- rnorm(3)
  y <- rnorm(3)

  expect_error(
    granger_causality_test(x = x, y = y, lag = 1),
    "Time series too short"
  )

  # With lag = 2, need at least 6 observations
  x <- rnorm(5)
  y <- rnorm(5)
  expect_error(
    granger_causality_test(x = x, y = y, lag = 2),
    "Time series too short"
  )
})

test_that("error on NA values in x", {
  x <- c(rnorm(50), NA, rnorm(49))
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = y),
    "contains NA"
  )
})

test_that("error on NA values in y", {
  x <- rnorm(100)
  y <- c(rnorm(50), NA, rnorm(49))

  expect_error(
    granger_causality_test(x = x, y = y),
    "contains NA"
  )
})

test_that("error on NaN values", {
  x <- c(rnorm(50), NaN, rnorm(49))
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = y),
    "contains NA, NaN"
  )
})

test_that("error on infinite values", {
  x <- c(rnorm(50), Inf, rnorm(49))
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = y),
    "infinite values"
  )

  x <- c(rnorm(50), -Inf, rnorm(49))
  expect_error(
    granger_causality_test(x = x, y = y),
    "infinite values"
  )
})

test_that("error on invalid lag values", {
  x <- rnorm(100)
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = y, lag = 0),
    "positive integer"
  )

  expect_error(
    granger_causality_test(x = x, y = y, lag = -1),
    "positive integer"
  )

  expect_error(
    granger_causality_test(x = x, y = y, lag = 1.5),
    "positive integer"
  )

  expect_error(
    granger_causality_test(x = x, y = y, lag = "one"),
    "positive integer"
  )

  expect_error(
    granger_causality_test(x = x, y = y, lag = c(1, 2)),
    "positive integer"
  )
})

test_that("error on invalid alpha values", {
  x <- rnorm(100)
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = y, alpha = 0),
    "between 0 and 1"
  )

  expect_error(
    granger_causality_test(x = x, y = y, alpha = 1),
    "between 0 and 1"
  )

  expect_error(
    granger_causality_test(x = x, y = y, alpha = -0.1),
    "between 0 and 1"
  )

  expect_error(
    granger_causality_test(x = x, y = y, alpha = 1.5),
    "between 0 and 1"
  )

  expect_error(
    granger_causality_test(x = x, y = y, alpha = "0.05"),
    "between 0 and 1"
  )
})

test_that("error on invalid test type", {
  x <- rnorm(100)
  y <- rnorm(100)

  expect_error(
    granger_causality_test(x = x, y = y, test = "Chi"),
    "only.*F.*supported"
  )

  expect_error(
    granger_causality_test(x = x, y = y, test = "t"),
    "only.*F.*supported"
  )
})

test_that("error on invalid .data argument", {
  expect_error(
    granger_causality_test(.data = "not a data frame", x = a, y = b),
    "must be a data frame"
  )

  expect_error(
    granger_causality_test(.data = list(a = 1:10, b = 11:20), x = a, y = b),
    "must be a data frame"
  )
})

test_that("error when column not found in data frame", {
  df <- data.frame(a = rnorm(100), b = rnorm(100))

  expect_error(
    granger_causality_test(df, x = a, y = nonexistent)
  )
})

test_that("works with minimum viable series length", {
  # For lag = 1, minimum is 2*1 + 2 = 4 observations
  set.seed(123)
  x <- rnorm(4)
  y <- rnorm(4)

  # Should work without error
  result <- granger_causality_test(x = x, y = y, lag = 1)
  expect_s3_class(result, "granger_result")
})
