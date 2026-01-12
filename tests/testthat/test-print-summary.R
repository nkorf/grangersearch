# Tests for print and summary methods

test_that("print method runs without error", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result <- granger_causality_test(x = x, y = y)

  expect_output(print(result), "Granger Causality Test")
  expect_output(print(result), "Observations:")
  expect_output(print(result), "Lag:")
})

test_that("print method returns object invisibly", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y)

  printed <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})

test_that("print method shows correct variable names", {
  skip_if_not_installed("tibble")

  df <- tibble::tibble(
    price = cumsum(rnorm(100)),
    volume = cumsum(rnorm(100))
  )

  result <- granger_causality_test(df, price, volume)

  expect_output(print(result), "price")
  expect_output(print(result), "volume")
})

test_that("print method shows p-values", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y)

  output <- capture.output(print(result))
  output_text <- paste(output, collapse = "\n")

  expect_match(output_text, "p = ")
})

test_that("summary method runs without error", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)

  result <- granger_causality_test(x = x, y = y)

  expect_output(summary(result), "Granger Causality Test - Detailed Summary")
  expect_output(summary(result), "Call:")
  expect_output(summary(result), "Results:")
  expect_output(summary(result), "Interpretation:")
})

test_that("summary method returns object invisibly", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y)

  printed <- capture.output(ret <- summary(result))
  expect_identical(ret, result)
})

test_that("summary shows bidirectional causality message when appropriate", {
  # Create data where both series cause each other
  set.seed(1)
  n <- 200
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- rnorm(1)
  y[1] <- rnorm(1)

  for (i in 2:n) {
    x[i] <- 0.7 * y[i-1] + rnorm(1, sd = 0.3)
    y[i] <- 0.7 * x[i-1] + rnorm(1, sd = 0.3)
  }

  result <- granger_causality_test(x = x, y = y)

  if (result$x_causes_y && result$y_causes_x) {
    expect_output(summary(result), "Bidirectional")
  }
})

test_that("summary shows unidirectional causality message", {
  set.seed(42)
  n <- 200
  x <- cumsum(rnorm(n))
  y <- c(0, 0.9 * x[1:(n-1)]) + rnorm(n, sd = 0.2)

  result <- granger_causality_test(x = x, y = y)

  if (result$x_causes_y && !result$y_causes_x) {
    expect_output(summary(result), "Unidirectional")
  }
})

test_that("summary shows test statistics", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y)

  expect_output(summary(result), "Test statistic:")
  expect_output(summary(result), "P-value:")
})

test_that("summary shows significance level used", {
  set.seed(123)
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))

  result <- granger_causality_test(x = x, y = y, alpha = 0.10)

  expect_output(summary(result), "0.100")
})
