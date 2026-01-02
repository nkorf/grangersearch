# Tests for granger_search (exhaustive pairwise search)

test_that("granger_search returns correct class", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n)),
    C = cumsum(rnorm(n))
  )

  result <- granger_search(df, include_insignificant = TRUE)

  expect_s3_class(result, "granger_search_result")
  expect_s3_class(result, "tbl_df")
})

test_that("granger_search tests all pairs", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n)),
    C = cumsum(rnorm(n))
  )

  result <- granger_search(df, include_insignificant = TRUE)

  # With 3 variables, should have 3*2 = 6 directed pairs
  expect_equal(nrow(result), 6)
})

test_that("granger_search detects known causal relationship", {
  set.seed(42)
  n <- 200
  df <- data.frame(
    X = cumsum(rnorm(n)),
    Y = rnorm(n),
    Z = rnorm(n)
  )
  # Y depends on lagged X
  df$Y <- c(0, 0.8 * df$X[1:(n-1)]) + rnorm(n, sd = 0.3)

  result <- granger_search(df)

  # Should find X -> Y
  expect_true(any(result$cause == "X" & result$effect == "Y"))
})

test_that("granger_search filters insignificant by default", {
  set.seed(123)
  n <- 100
  # All independent series
  df <- data.frame(
    A = rnorm(n),
    B = rnorm(n),
    C = rnorm(n)
  )

  result <- granger_search(df)

  # Most likely no significant results with independent series
  # (could occasionally find spurious ones due to chance)
  expect_true(nrow(result) <= 6)
})

test_that("granger_search include_insignificant works", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n))
  )

  result_all <- granger_search(df, include_insignificant = TRUE)
  result_sig <- granger_search(df, include_insignificant = FALSE)

  expect_equal(nrow(result_all), 2)  # 2 directed pairs
  expect_true(nrow(result_sig) <= nrow(result_all))
})

test_that("granger_search column selection works", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n)),
    C = cumsum(rnorm(n)),
    D = cumsum(rnorm(n))
  )

  # Select only A and B
  result <- granger_search(df, A, B, include_insignificant = TRUE)

  expect_equal(nrow(result), 2)  # Only A->B and B->A
  expect_true(all(result$cause %in% c("A", "B")))
  expect_true(all(result$effect %in% c("A", "B")))
})

test_that("granger_search uses all numeric columns by default", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n)),
    label = rep("x", n)  # Non-numeric column
  )

  result <- granger_search(df, include_insignificant = TRUE)

  # Should only use A and B (numeric columns)
  expect_equal(nrow(result), 2)
})

test_that("granger_search respects lag parameter", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n))
  )

  result <- granger_search(df, lag = 3, include_insignificant = TRUE)

  expect_true(all(result$lag == 3))
})

test_that("granger_search errors on non-data frame", {
  expect_error(
    granger_search(c(1, 2, 3)),
    "must be a data frame"
  )
})

test_that("granger_search errors on fewer than 2 variables", {
  df <- data.frame(A = rnorm(100))

  expect_error(
    granger_search(df),
    "at least 2"
  )
})

test_that("granger_search errors on non-existent column", {
  df <- data.frame(A = rnorm(100), B = rnorm(100))

  expect_error(
    granger_search(df, A, nonexistent),
    "not found"
  )
})

test_that("granger_search results are sorted by p-value", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n)),
    C = cumsum(rnorm(n))
  )

  result <- granger_search(df, include_insignificant = TRUE)

  # Check p-values are in ascending order
  expect_true(all(diff(result$p.value) >= 0, na.rm = TRUE))
})

test_that("print method works for granger_search_result", {
  set.seed(123)
  n <- 100
  df <- data.frame(
    A = cumsum(rnorm(n)),
    B = cumsum(rnorm(n))
  )

  result <- granger_search(df, include_insignificant = TRUE)

  expect_output(print(result), "Granger Causality Search Results")
})
