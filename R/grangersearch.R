#' Perform a Granger Causality Test on Two Time Series
#'
#' Tests whether one time series Granger-causes another and vice versa.
#' A variable X is said to Granger-cause Y if past values of X help predict Y
#' beyond what past values of Y alone provide.
#'
#' @param .data A data frame, tibble, or NULL. If provided, `x` and `y` are
#'   evaluated as column names within this data frame (tidyverse-style).
#' @param x Either a numeric vector/time series, or (if `.data` is provided)
#'   an unquoted column name.
#' @param y Either a numeric vector/time series of the same length as `x`, or
#'   (if `.data` is provided) an unquoted column name.
#' @param lag Integer. The lag order for the VAR model. Default is 1.
#' @param alpha Numeric. Significance level for the causality test (between 0 and 1).
#'   Default is 0.05.
#' @param test Character. Type of test to perform. Currently only "F" (F-test) is
#'   supported. Default is "F".
#'
#' @return An object of class `granger_result` containing:
#' \describe{
#'   \item{x_causes_y}{Logical. TRUE if X Granger-causes Y at the specified alpha level.}
#'   \item{y_causes_x}{Logical. TRUE if Y Granger-causes X at the specified alpha level.}
#'   \item{p_value_xy}{Numeric. P-value for the test of X causing Y.}
#'   \item{p_value_yx}{Numeric. P-value for the test of Y causing X.}
#'   \item{test_statistic_xy}{Numeric. Test statistic for X causing Y.}
#'   \item{test_statistic_yx}{Numeric. Test statistic for Y causing X.}
#'   \item{lag}{Integer. The lag order used.}
#'   \item{alpha}{Numeric. The significance level used.}
#'   \item{test}{Character. The test type used.}
#'   \item{n}{Integer. Number of observations.}
#'   \item{x_name}{Character. Name of the X variable.}
#'   \item{y_name}{Character. Name of the Y variable.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The Granger causality test is based on the idea that if X causes Y, then past
#' values of X should contain information that helps predict Y above and beyond
#' the information contained in past values of Y alone.
#'
#' This function fits Vector Autoregressive (VAR) models using the \pkg{vars}
#' package and performs F-tests to compare restricted and unrestricted models.
#' The test is performed in both directions to detect unidirectional or
#' bidirectional causality.
#'
#' Note that Granger causality is a statistical concept based on prediction and
#' temporal precedence. It does not necessarily imply true causal mechanisms.
#'
#' @section Tidyverse Compatibility:
#' This function supports tidyverse-style syntax:
#' \itemize{
#'   \item Pipe-friendly: use with `%>%` or `|>`
#'   \item NSE column selection: pass unquoted column names when using a data frame
#'   \item Use [tidy.granger_result()] to get a tibble of results
#'   \item Use [glance.granger_result()] for model-level summary
#' }
#'
#' @examples
#' # Vector-based usage
#' set.seed(123)
#' n <- 100
#' x <- cumsum(rnorm(n))
#' y <- c(0, x[1:(n-1)]) + rnorm(n, sd = 0.5)
#'
#' result <- granger_causality_test(x = x, y = y)
#' print(result)
#'
#' # Tidyverse-style with data frame
#' library(tibble)
#' df <- tibble(
#'   price = cumsum(rnorm(100)),
#'   volume = c(0, cumsum(rnorm(99)))
#' )
#'
#' # Using pipe and column names
#' df |> granger_causality_test(price, volume)
#'
#' # Get tidy results as tibble
#' result |> tidy()
#'
#' # Different lag order
#' df |> granger_causality_test(price, volume, lag = 2)
#'
#' @references
#' Granger, C. W. J. (1969). Investigating Causal Relations by Econometric Models
#' and Cross-spectral Methods. \emph{Econometrica}, 37(3), 424-438.
#'
#' @seealso
#' \code{\link[vars]{VAR}} for the underlying VAR model,
#' \code{\link[vars]{causality}} for an alternative implementation,
#' [tidy.granger_result()] for tidying results.
#'
#' @importFrom rlang enquo eval_tidy as_label
#' @importFrom stats complete.cases
#' @export
granger_causality_test <- function(.data = NULL, x, y, lag = 1, alpha = 0.05, test = "F") {
  call <- match.call()

  # Handle tidyverse-style NSE
  x_quo <- rlang::enquo(x)
  y_quo <- rlang::enquo(y)

  if (!is.null(.data)) {
    # Data frame provided - evaluate column names
    if (!is.data.frame(.data)) {
      stop("`.data` must be a data frame or tibble.", call. = FALSE)
    }
    x_vec <- rlang::eval_tidy(x_quo, .data)
    y_vec <- rlang::eval_tidy(y_quo, .data)
    x_name <- rlang::as_label(x_quo)
    y_name <- rlang::as_label(y_quo)
  } else {
    # Direct vectors provided
    x_vec <- rlang::eval_tidy(x_quo)
    y_vec <- rlang::eval_tidy(y_quo)
    x_name <- rlang::as_label(x_quo)
    y_name <- rlang::as_label(y_quo)
  }

  # Input validation
  if (!is.numeric(x_vec)) {
    stop(sprintf("`%s` must be numeric.", x_name), call. = FALSE)
  }
  if (!is.numeric(y_vec)) {
    stop(sprintf("`%s` must be numeric.", y_name), call. = FALSE)
  }
  if (length(x_vec) != length(y_vec)) {
    stop(sprintf("`%s` and `%s` must have the same length.", x_name, y_name), call. = FALSE)
  }

  # Validate parameters before using them
  if (!is.numeric(lag) || length(lag) != 1 || lag < 1 || lag != as.integer(lag)) {
    stop("`lag` must be a positive integer.", call. = FALSE)
  }
  lag <- as.integer(lag)

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a number between 0 and 1 (exclusive).", call. = FALSE)
  }

  if (!identical(test, "F")) {
    stop("Currently only `test = 'F'` is supported.", call. = FALSE)
  }

  n <- length(x_vec)
  min_length <- 2 * lag + 2
  if (n < min_length) {
    stop(sprintf("Time series too short. Need at least %d observations for lag = %d.",
                 min_length, lag), call. = FALSE)
  }

  if (any(is.na(x_vec)) || any(is.nan(x_vec)) || any(is.infinite(x_vec))) {
    stop(sprintf("`%s` contains NA, NaN, or infinite values.", x_name), call. = FALSE)
  }
  if (any(is.na(y_vec)) || any(is.nan(y_vec)) || any(is.infinite(y_vec))) {
    stop(sprintf("`%s` contains NA, NaN, or infinite values.", y_name), call. = FALSE)
  }

  # Combine into matrix for VAR modeling
  data_matrix <- cbind(X = x_vec, Y = y_vec)

  # Fit VAR model and perform causality tests
  tryCatch({
    var_model <- vars::VAR(data_matrix, p = lag, type = "const")

    # Test if X Granger-causes Y
    causality_xy <- vars::causality(var_model, cause = "X")
    p_value_xy <- causality_xy$Granger$p.value
    test_stat_xy <- causality_xy$Granger$statistic

    # Test if Y Granger-causes X
    causality_yx <- vars::causality(var_model, cause = "Y")
    p_value_yx <- causality_yx$Granger$p.value
    test_stat_yx <- causality_yx$Granger$statistic
  }, error = function(e) {
    stop(sprintf("VAR model fitting failed: %s", e$message), call. = FALSE)
  })

  # Build result object
  result <- list(
    x_causes_y = p_value_xy < alpha,
    y_causes_x = p_value_yx < alpha,
    p_value_xy = as.numeric(p_value_xy),
    p_value_yx = as.numeric(p_value_yx),
    test_statistic_xy = as.numeric(test_stat_xy),
    test_statistic_yx = as.numeric(test_stat_yx),
    lag = lag,
    alpha = alpha,
    test = test,
    n = n,
    x_name = x_name,
    y_name = y_name,
    call = call
  )

  class(result) <- "granger_result"
  result
}


#' Print Method for granger_result Objects
#'
#' @param x A `granger_result` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.granger_result <- function(x, ...) {
  cat("\nGranger Causality Test\n")
  cat("======================\n\n")

  cat(sprintf("Observations: %d, Lag order: %d, Significance level: %.3f\n\n",
              x$n, x$lag, x$alpha))

  # X -> Y result
  xy_result <- if (x$x_causes_y) {
    sprintf("%s Granger-causes %s", x$x_name, x$y_name)
  } else {
    sprintf("%s does not Granger-cause %s", x$x_name, x$y_name)
  }
  cat(sprintf("%s -> %s: %s (p = %.4f)\n", x$x_name, x$y_name, xy_result, x$p_value_xy))

  # Y -> X result
  yx_result <- if (x$y_causes_x) {
    sprintf("%s Granger-causes %s", x$y_name, x$x_name)
  } else {
    sprintf("%s does not Granger-cause %s", x$y_name, x$x_name)
  }
  cat(sprintf("%s -> %s: %s (p = %.4f)\n", x$y_name, x$x_name, yx_result, x$p_value_yx))

  cat("\n")
  invisible(x)
}


#' Summary Method for granger_result Objects
#'
#' @param object A `granger_result` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the object.
#'
#' @export
summary.granger_result <- function(object, ...) {
  cat("\nGranger Causality Test - Detailed Summary\n")
  cat("==========================================\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat(sprintf("Variables: %s, %s\n", object$x_name, object$y_name))
  cat(sprintf("Number of observations: %d\n", object$n))
  cat(sprintf("VAR lag order: %d\n", object$lag))
  cat(sprintf("Test type: %s-test\n", object$test))
  cat(sprintf("Significance level (alpha): %.3f\n\n", object$alpha))

  cat("Results:\n")
  cat("--------\n\n")

  cat(sprintf("Test: %s Granger-causes %s\n", object$x_name, object$y_name))
  cat(sprintf("  Test statistic: %.4f\n", object$test_statistic_xy))
  cat(sprintf("  P-value: %.6f\n", object$p_value_xy))
  cat(sprintf("  Conclusion: %s\n\n",
              if (object$x_causes_y) {
                sprintf("REJECT null (%s causes %s)", object$x_name, object$y_name)
              } else {
                sprintf("FAIL TO REJECT null (%s does not cause %s)", object$x_name, object$y_name)
              }))

  cat(sprintf("Test: %s Granger-causes %s\n", object$y_name, object$x_name))
  cat(sprintf("  Test statistic: %.4f\n", object$test_statistic_yx))
  cat(sprintf("  P-value: %.6f\n", object$p_value_yx))
  cat(sprintf("  Conclusion: %s\n\n",
              if (object$y_causes_x) {
                sprintf("REJECT null (%s causes %s)", object$y_name, object$x_name)
              } else {
                sprintf("FAIL TO REJECT null (%s does not cause %s)", object$y_name, object$x_name)
              }))

  # Overall interpretation
  cat("Interpretation:\n")
  cat("---------------\n")
  if (object$x_causes_y && object$y_causes_x) {
    cat(sprintf("Bidirectional causality detected: %s and %s Granger-cause each other.\n",
                object$x_name, object$y_name))
  } else if (object$x_causes_y) {
    cat(sprintf("Unidirectional causality: %s Granger-causes %s (but not vice versa).\n",
                object$x_name, object$y_name))
  } else if (object$y_causes_x) {
    cat(sprintf("Unidirectional causality: %s Granger-causes %s (but not vice versa).\n",
                object$y_name, object$x_name))
  } else {
    cat("No Granger causality detected in either direction.\n")
  }
  cat("\n")

  invisible(object)
}


#' Tidy a granger_result Object
#'
#' Returns a tibble with one row per direction tested, containing test results.
#' Compatible with the broom package conventions.
#'
#' @param x A `granger_result` object.
#' @param ... Additional arguments (ignored).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{direction}{Character. The causal direction tested (e.g., "x -> y").}
#'   \item{cause}{Character. The name of the potential cause variable.}
#'   \item{effect}{Character. The name of the potential effect variable.}
#'   \item{statistic}{Numeric. The F-test statistic.}
#'   \item{p.value}{Numeric. The p-value of the test.}
#'   \item{significant}{Logical. Whether the result is significant at the alpha level.}
#' }
#'
#' @examples
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#' y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)
#' result <- granger_causality_test(x = x, y = y)
#' tidy(result)
#'
#' @importFrom tibble tibble
#' @method tidy granger_result
#' @export
tidy.granger_result <- function(x, ...) {
tibble::tibble(
    direction = c(
      sprintf("%s -> %s", x$x_name, x$y_name),
      sprintf("%s -> %s", x$y_name, x$x_name)
    ),
    cause = c(x$x_name, x$y_name),
    effect = c(x$y_name, x$x_name),
    statistic = c(x$test_statistic_xy, x$test_statistic_yx),
    p.value = c(x$p_value_xy, x$p_value_yx),
    significant = c(x$x_causes_y, x$y_causes_x)
  )
}


#' Glance at a granger_result Object
#'
#' Returns a tibble with a single row containing model-level summary statistics.
#' Compatible with the broom package conventions.
#'
#' @param x A `granger_result` object.
#' @param ... Additional arguments (ignored).
#'
#' @return A tibble with one row and columns:
#' \describe{
#'   \item{nobs}{Integer. Number of observations.}
#'   \item{lag}{Integer. VAR lag order used.}
#'   \item{alpha}{Numeric. Significance level used.}
#'   \item{test}{Character. Test type used.}
#'   \item{bidirectional}{Logical. TRUE if causality detected in both directions.}
#'   \item{x_causes_y}{Logical. TRUE if x Granger-causes y.}
#'   \item{y_causes_x}{Logical. TRUE if y Granger-causes x.}
#' }
#'
#' @examples
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#' y <- c(0, x[1:99]) + rnorm(100, sd = 0.5)
#' result <- granger_causality_test(x = x, y = y)
#' glance(result)
#'
#' @importFrom tibble tibble
#' @method glance granger_result
#' @export
glance.granger_result <- function(x, ...) {
  tibble::tibble(
    nobs = x$n,
    lag = x$lag,
    alpha = x$alpha,
    test = x$test,
    bidirectional = x$x_causes_y && x$y_causes_x,
    x_causes_y = x$x_causes_y,
    y_causes_x = x$y_causes_x
  )
}


#' Exhaustive Pairwise Granger Causality Search
#'
#' Performs Granger causality tests on all pairwise combinations of variables
#' in a dataset. This is the core "search" functionality of the package,
#' enabling discovery of causal relationships among multiple time series.
#'
#' @param .data A data frame or tibble containing the time series variables.
#' @param ... <[`tidy-select`][dplyr::dplyr_tidy_select]> Columns to include
#'   in the analysis. If empty, all numeric columns are used.
#' @param lag Integer or integer vector. The lag order(s) for VAR models.
#'   If a vector (e.g., `1:4`), tests are performed at each lag and the
#'   best (lowest p-value) result is returned for each pair. Default is 1.
#' @param alpha Numeric. Significance level for hypothesis testing. Default is 0.05.
#' @param test Character. Test type, currently only "F" supported. Default is "F".
#' @param include_insignificant Logical. If FALSE (default), only return
#'   significant causal relationships. If TRUE, return all pairwise results.
#'
#' @return A tibble with one row per directed pair tested, containing:
#' \describe{
#'   \item{cause}{Character. The potential cause variable name.}
#'   \item{effect}{Character. The potential effect variable name.}
#'   \item{statistic}{Numeric. The F-test statistic.}
#'   \item{p.value}{Numeric. The p-value of the test.}
#'   \item{significant}{Logical. Whether the result is significant at alpha.}
#'   \item{lag}{Integer. The lag order used (best lag if multiple were tested).}
#' }
#'
#' @details
#' This function tests all \eqn{n(n-1)} directed pairs for \eqn{n} variables.
#' For each pair (X, Y), it tests whether X Granger-causes Y.
#'
#' When multiple lags are specified (e.g., `lag = 1:4`), the function tests
#' each pair at every lag and returns the result with the lowest p-value.
#' This is useful for discovering the optimal lag structure.
#'
#' The function is useful for exploratory analysis when you have multiple
#' time series and want to discover which variables have predictive relationships.
#'
#' @section Multiple Testing:
#' When testing many pairs (and especially many lags), consider adjusting for
#' multiple comparisons. The returned p-values are unadjusted. You can apply
#' corrections such as Bonferroni or Benjamini-Hochberg using [stats::p.adjust()].
#'
#' @examples
#' # Create dataset with multiple time series
#' set.seed(123)
#' n <- 100
#' df <- data.frame(
#'   A = cumsum(rnorm(n)),
#'   B = cumsum(rnorm(n)),
#'   C = cumsum(rnorm(n))
#' )
#' # B is caused by lagged A
#' df$B <- c(0, 0.7 * df$A[1:(n-1)]) + rnorm(n, sd = 0.5)
#'
#' # Search for all causal relationships
#' granger_search(df)
#'
#' # Include all results, not just significant ones
#' granger_search(df, include_insignificant = TRUE)
#'
#' # Select specific columns
#' granger_search(df, A, B)
#'
#' # Search across multiple lags (returns best lag for each pair)
#' granger_search(df, lag = 1:4)
#'
#' # Search with specific lag
#' granger_search(df, lag = 2)
#'
#' @seealso [granger_causality_test()] for testing a single pair.
#'
#' @importFrom tibble tibble
#' @importFrom rlang enquos eval_tidy as_label
#' @export
granger_search <- function(.data, ..., lag = 1, alpha = 0.05, test = "F",
                           include_insignificant = FALSE) {

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  # Validate lag parameter
  if (!is.numeric(lag) || any(lag < 1) || any(lag != as.integer(lag))) {
    stop("`lag` must be positive integer(s).", call. = FALSE)
  }
  lag <- as.integer(lag)

  # Handle column selection
  cols <- rlang::enquos(...)

  if (length(cols) == 0) {
    # No columns specified - use all numeric columns
    var_names <- names(.data)[sapply(.data, is.numeric)]
    if (length(var_names) < 2) {
      stop("Need at least 2 numeric columns for pairwise testing.", call. = FALSE)
    }
  } else {
    # Evaluate tidyselect
    var_names <- character()
    for (col in cols) {
      var_names <- c(var_names, rlang::as_label(col))
    }
    # Verify columns exist and are numeric
    for (v in var_names) {
      if (!v %in% names(.data)) {
        stop(sprintf("Column `%s` not found in data.", v), call. = FALSE)
      }
      if (!is.numeric(.data[[v]])) {
        stop(sprintf("Column `%s` must be numeric.", v), call. = FALSE)
      }
    }
  }

  n_vars <- length(var_names)
  if (n_vars < 2) {
    stop("Need at least 2 variables for pairwise testing.", call. = FALSE)
  }

  # Generate all ordered pairs
  pairs <- expand.grid(cause = var_names, effect = var_names,
                       stringsAsFactors = FALSE)
  pairs <- pairs[pairs$cause != pairs$effect, ]

  # Test each pair (across all lags if multiple specified)
  results <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    cause_var <- pairs$cause[i]
    effect_var <- pairs$effect[i]

    # Get vectors
    x_vec <- .data[[cause_var]]
    y_vec <- .data[[effect_var]]

    # Test at each lag and keep best result
    best_result <- NULL
    best_p <- Inf

    for (p in lag) {
      tryCatch({
        result <- granger_causality_test(
          x = x_vec, y = y_vec,
          lag = p, alpha = alpha, test = test
        )

        if (result$p_value_xy < best_p) {
          best_p <- result$p_value_xy
          best_result <- tibble::tibble(
            cause = cause_var,
            effect = effect_var,
            statistic = result$test_statistic_xy,
            p.value = result$p_value_xy,
            significant = result$x_causes_y,
            lag = p
          )
        }
      }, error = function(e) {
        # Skip this lag if it fails
      })
    }

    if (is.null(best_result)) {
      # All lags failed
      results[[i]] <- tibble::tibble(
        cause = cause_var,
        effect = effect_var,
        statistic = NA_real_,
        p.value = NA_real_,
        significant = NA,
        lag = lag[1]
      )
      warning(sprintf("All tests failed for %s -> %s", cause_var, effect_var),
              call. = FALSE)
    } else {
      results[[i]] <- best_result
    }
  }

  # Combine results
  output <- do.call(rbind, results)

  # Filter if requested
  if (!include_insignificant) {
    output <- output[output$significant == TRUE & !is.na(output$significant), ]
  }

  # Sort by p-value
  output <- output[order(output$p.value), ]
  rownames(output) <- NULL

  # Store metadata for potential use
  attr(output, "lags_tested") <- lag
  attr(output, "alpha") <- alpha


  # Add class for S3 methods - must be first for proper dispatch
  class(output) <- c("granger_search_result", class(output))

  output
}


#' Lag Selection Analysis for Granger Causality
#'
#' Analyzes how Granger causality test results change across different lag orders.
#' Returns detailed results for all lag-pair combinations, useful for optimal
#' lag selection and visualization.
#'
#' @param .data A data frame or tibble containing the time series variables.
#' @param ... <[`tidy-select`][dplyr::dplyr_tidy_select]> Columns to include.
#'   If empty, all numeric columns are used.
#' @param lag Integer vector. The lag orders to test. Default is `1:4`.
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param test Character. Test type. Default is "F".
#'
#' @return A tibble with one row per (cause, effect, lag) combination:
#' \describe{
#'   \item{cause}{Character. The potential cause variable.}
#'   \item{effect}{Character. The potential effect variable.}
#'   \item{lag}{Integer. The lag order tested.}
#'   \item{statistic}{Numeric. The F-test statistic.}
#'   \item{p.value}{Numeric. The p-value.}
#'   \item{significant}{Logical. Whether significant at alpha.}
#' }
#'
#' @details
#' Unlike [granger_search()] which returns only the best lag for each pair,
#' this function returns results for all lag values tested. This is useful for:
#' \itemize{
#'   \item Visualizing how p-values change with lag order
#'   \item Selecting the optimal lag for each relationship
#'   \item Understanding the temporal dynamics of causality
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' df <- data.frame(
#'   A = cumsum(rnorm(n)),
#'   B = cumsum(rnorm(n))
#' )
#' df$B <- c(0, 0.7 * df$A[1:(n-1)]) + rnorm(n, sd = 0.5)
#'
#' # Get results for lags 1 through 5
#' lag_results <- granger_lag_select(df, lag = 1:5)
#'
#' # Can be used with ggplot2 for visualization
#' # library(ggplot2)
#' # ggplot(lag_results, aes(x = lag, y = p.value, color = paste(cause, "->", effect))) +
#' #   geom_line() + geom_point() +
#' #   geom_hline(yintercept = 0.05, linetype = "dashed") +
#' #   labs(title = "P-values by Lag Order", color = "Direction")
#'
#' @seealso [granger_search()] for getting best results across lags,
#'   [plot.granger_lag_select()] for built-in visualization.
#'
#' @importFrom tibble tibble
#' @export
granger_lag_select <- function(.data, ..., lag = 1:4, alpha = 0.05, test = "F") {

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(lag) || any(lag < 1) || any(lag != as.integer(lag))) {
    stop("`lag` must be positive integer(s).", call. = FALSE)
  }
  lag <- as.integer(lag)

  # Handle column selection
  cols <- rlang::enquos(...)

  if (length(cols) == 0) {
    var_names <- names(.data)[sapply(.data, is.numeric)]
    if (length(var_names) < 2) {
      stop("Need at least 2 numeric columns.", call. = FALSE)
    }
  } else {
    var_names <- character()
    for (col in cols) {
      var_names <- c(var_names, rlang::as_label(col))
    }
    for (v in var_names) {
      if (!v %in% names(.data)) {
        stop(sprintf("Column `%s` not found.", v), call. = FALSE)
      }
      if (!is.numeric(.data[[v]])) {
        stop(sprintf("Column `%s` must be numeric.", v), call. = FALSE)
      }
    }
  }

  # Generate all pairs
  pairs <- expand.grid(cause = var_names, effect = var_names,
                       stringsAsFactors = FALSE)
  pairs <- pairs[pairs$cause != pairs$effect, ]

  # Test each pair at each lag
  results <- list()

  for (i in seq_len(nrow(pairs))) {
    cause_var <- pairs$cause[i]
    effect_var <- pairs$effect[i]
    x_vec <- .data[[cause_var]]
    y_vec <- .data[[effect_var]]

    for (p in lag) {
      tryCatch({
        result <- granger_causality_test(
          x = x_vec, y = y_vec,
          lag = p, alpha = alpha, test = test
        )

        results[[length(results) + 1]] <- tibble::tibble(
          cause = cause_var,
          effect = effect_var,
          lag = p,
          statistic = result$test_statistic_xy,
          p.value = result$p_value_xy,
          significant = result$x_causes_y
        )
      }, error = function(e) {
        results[[length(results) + 1]] <<- tibble::tibble(
          cause = cause_var,
          effect = effect_var,
          lag = p,
          statistic = NA_real_,
          p.value = NA_real_,
          significant = NA
        )
      })
    }
  }

  output <- do.call(rbind, results)
  attr(output, "alpha") <- alpha
  class(output) <- c("granger_lag_select", class(output))
  output
}


#' Plot Lag Selection Results
#'
#' Creates a visualization of p-values across different lag orders for
#' Granger causality tests.
#'
#' @param x A `granger_lag_select` object from [granger_lag_select()].
#' @param ... Additional arguments (ignored).
#'
#' @return A base R plot (invisibly returns the input).
#'
#' @details
#' This function creates a line plot showing how p-values change across
#' different lag orders for each directed pair. A horizontal dashed line
#' indicates the significance threshold (alpha).
#'
#' For more customized plots, use the data directly with ggplot2.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(A = cumsum(rnorm(100)), B = cumsum(rnorm(100)))
#' df$B <- c(0, 0.7 * df$A[1:99]) + rnorm(100, sd = 0.5)
#'
#' lag_results <- granger_lag_select(df, lag = 1:5)
#' plot(lag_results)
#'
#' @export
plot.granger_lag_select <- function(x, ...) {
  alpha <- attr(x, "alpha") %||% 0.05

  # Get unique directions
  x$direction <- paste(x$cause, "->", x$effect)
  directions <- unique(x$direction)
  n_dir <- length(directions)

  # Set up colors
  colors <- grDevices::rainbow(n_dir)

  # Create plot
  graphics::plot(
    range(x$lag), c(0, max(x$p.value, na.rm = TRUE) * 1.1),
    type = "n",
    xlab = "Lag Order",
    ylab = "P-value",
    main = "Granger Causality: P-values by Lag Order"
  )

  # Add horizontal line at alpha
  graphics::abline(h = alpha, lty = 2, col = "gray50")

  # Plot each direction
  for (i in seq_along(directions)) {
    d <- directions[i]
    subset_data <- x[x$direction == d, ]
    graphics::lines(subset_data$lag, subset_data$p.value, col = colors[i], lwd = 2)
    graphics::points(subset_data$lag, subset_data$p.value, col = colors[i], pch = 19)
  }

  # Add legend
  graphics::legend(
    "topright",
    legend = directions,
    col = colors,
    lwd = 2,
    pch = 19,
    cex = 0.8,
    bg = "white"
  )

  invisible(x)
}

# Helper for NULL default
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Print Method for granger_search_result Objects
#'
#' @param x A `granger_search_result` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.granger_search_result <- function(x, ...) {
  n_sig <- sum(x$significant, na.rm = TRUE)
  n_total <- nrow(x)

  cat("\nGranger Causality Search Results\n")
  cat("=================================\n\n")

  if (n_total == 0) {
    cat("No significant causal relationships found.\n")
  } else {
    cat(sprintf("Found %d significant relationship(s):\n\n", n_sig))
    # Print as tibble
    print(tibble::as_tibble(x))
  }


  cat("\n")
  invisible(x)
}


#' Plot Granger Causality Matrix
#'
#' Creates a heatmap-style matrix visualization of Granger causality relationships.
#' Displays two panels showing both directions of causality testing.
#'
#' @param x A `granger_search_result` object from [granger_search()].
#' @param type Character. Type of values to display:
#'   \itemize{
#'     \item `"pvalue"` (default): Show -log10(p-value), with higher values indicating
#'       stronger evidence of causality.
#'     \item `"significance"`: Show binary significant/not significant.
#'     \item `"statistic"`: Show the F-test statistic values.
#'   }
#' @param show_values Logical. If TRUE, display numeric values in cells. Default is TRUE.
#' @param gradient Logical. If TRUE (default), use gradient coloring based on values.
#'   If FALSE, cells are colored simply as significant (colored) or not significant (gray),
#'   and actual p-values are displayed in cells when `type = "pvalue"`.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @details
#' The visualization shows two panels:
#' \itemize{
#'   \item Left panel: Tests whether row variable Granger-causes column variable
#'   \item Right panel: Tests whether column variable Granger-causes row variable
#' }
#'
#' For the `"pvalue"` type with `gradient = TRUE`, values are shown as -log10(p-value),
#' so larger values indicate stronger evidence of Granger causality. When `gradient = FALSE`,
#' actual p-values are displayed in cells and coloring is binary (significant vs not).
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' df <- data.frame(
#'   gdp = cumsum(rnorm(n)),
#'   consumption = cumsum(rnorm(n)),
#'   investment = cumsum(rnorm(n)),
#'   employment = cumsum(rnorm(n))
#' )
#' # Add some causal structure
#' df$consumption <- c(0, 0.5 * df$gdp[1:(n-1)]) + rnorm(n, sd = 0.5)
#' df$employment <- c(0, 0.3 * df$gdp[1:(n-1)]) + rnorm(n, sd = 0.5)
#'
#' # Run exhaustive search (include all results for complete matrix)
#' results <- granger_search(df, include_insignificant = TRUE)
#'
#' # Plot as matrix with gradient
#' plot(results)
#'
#' # Plot without gradient, showing p-values
#' plot(results, gradient = FALSE)
#'
#' # Show binary significance
#' plot(results, type = "significance")
#'
#' @seealso [granger_search()] for running the exhaustive search.
#'
#' @export
plot.granger_search_result <- function(x, type = c("pvalue", "significance", "statistic"),
                                        show_values = TRUE, gradient = TRUE, ...) {
  type <- match.arg(type)
  alpha <- attr(x, "alpha") %||% 0.05

  # Get unique variables
  all_vars <- unique(c(x$cause, x$effect))
  n_vars <- length(all_vars)

  # Create matrices for X -> Y
  # mat_xy holds display values (for coloring)
  # mat_pval holds actual p-values (for cell labels when gradient=FALSE)
  mat_xy <- matrix(NA, nrow = n_vars, ncol = n_vars,
                   dimnames = list(all_vars, all_vars))
  mat_pval <- matrix(NA, nrow = n_vars, ncol = n_vars,
                     dimnames = list(all_vars, all_vars))
  mat_sig <- matrix(NA, nrow = n_vars, ncol = n_vars,
                    dimnames = list(all_vars, all_vars))

  # Fill matrices
  for (i in seq_len(nrow(x))) {
    cause <- x$cause[i]
    effect <- x$effect[i]

    mat_pval[cause, effect] <- x$p.value[i]
    mat_sig[cause, effect] <- as.numeric(x$significant[i])

    if (type == "pvalue") {
      pval <- x$p.value[i]
      if (!is.na(pval) && pval > 0) {
        mat_xy[cause, effect] <- -log10(pval)
      }
    } else if (type == "significance") {
      mat_xy[cause, effect] <- as.numeric(x$significant[i])
    } else {
      mat_xy[cause, effect] <- x$statistic[i]
    }
  }

  # Create transposed matrices for Y -> X view
  mat_yx <- t(mat_xy)
  mat_pval_yx <- t(mat_pval)
  mat_sig_yx <- t(mat_sig)

  # Calculate max value for color scaling (same scale for both)
  max_val <- max(mat_xy, na.rm = TRUE)
  if (is.na(max_val) || max_val <= 0) max_val <- 1

  # Set up layout: two panels (+ optional color legend)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (gradient) {
    graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4.5, 3.5, 1.2))
  } else {
    # No legend needed when gradient is off
    # Left panel wider to accommodate y-axis labels
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.3, 1))
  }

  # Color palette
  if (!gradient || type == "significance") {
    colors <- c("gray90", "steelblue3")
    breaks <- c(-0.5, 0.5, 1.5)
    n_colors <- 2
  } else {
    n_colors <- 100
    colors <- grDevices::colorRampPalette(c("gray95", "lightblue", "steelblue", "darkblue"))(n_colors)
    breaks <- seq(0, max_val * 1.01, length.out = n_colors + 1)
  }

  # Helper function to draw a single heatmap panel
  draw_panel <- function(mat, mat_pval_panel, mat_sig_panel, title, show_y_labels = TRUE) {
    mat_plot <- mat[rev(seq_len(nrow(mat))), , drop = FALSE]
    mat_pval_plot <- mat_pval_panel[rev(seq_len(nrow(mat_pval_panel))), , drop = FALSE]
    mat_sig_plot <- mat_sig_panel[rev(seq_len(nrow(mat_sig_panel))), , drop = FALSE]

    # For non-gradient mode, use significance for coloring
    if (!gradient && type == "pvalue") {
      mat_color <- mat_sig_plot
    } else {
      mat_color <- mat_plot
    }

    graphics::image(
      x = seq_len(ncol(mat_color)),
      y = seq_len(nrow(mat_color)),
      z = t(mat_color),
      col = colors,
      breaks = breaks,
      axes = FALSE,
      xlab = "",
      ylab = ""
    )

    # Title
    graphics::title(main = title, cex.main = 1.1, font.main = 2, line = 0.5)

    # X-axis labels (angled) - positioned closer to axis
    graphics::axis(1, at = seq_len(n_vars), labels = FALSE, tick = TRUE)
    graphics::text(x = seq_len(n_vars), y = graphics::par("usr")[3] - 0.15,
                   labels = colnames(mat_plot),
                   srt = 45, adj = 1, xpd = TRUE, cex = 0.85)

    # Y-axis labels (only if requested)
    if (show_y_labels) {
      graphics::axis(2, at = seq_len(n_vars), labels = rownames(mat_plot),
                     las = 1, cex.axis = 0.85, tick = TRUE)
    } else {
      graphics::axis(2, at = seq_len(n_vars), labels = FALSE, tick = TRUE)
    }

    # Grid lines
    graphics::abline(h = seq(0.5, n_vars + 0.5, 1), col = "white", lwd = 1.5)
    graphics::abline(v = seq(0.5, n_vars + 0.5, 1), col = "white", lwd = 1.5)

    # Values in cells
    if (show_values && n_vars <= 8) {
      for (i in seq_len(n_vars)) {
        for (j in seq_len(n_vars)) {
          sig_val <- mat_sig_plot[i, j]
          pval <- mat_pval_plot[i, j]
          val <- mat_plot[i, j]

          if (!is.na(pval)) {
            # Determine text color based on background
            if (!gradient || type == "significance") {
              text_col <- if (!is.na(sig_val) && sig_val > 0.5) "white" else "gray30"
            } else {
              text_col <- if (!is.na(val) && val > max_val * 0.4) "white" else "gray30"
            }

            # Determine label
            if (type == "significance") {
              label <- if (!is.na(sig_val) && sig_val > 0.5) "\u2713" else ""
            } else if (!gradient && type == "pvalue") {
              # Show actual p-value when gradient is off
              label <- sprintf("%.3f", pval)
            } else if (type == "pvalue") {
              label <- sprintf("%.1f", val)
            } else {
              label <- sprintf("%.1f", val)
            }

            graphics::text(j, i, label, col = text_col, cex = 0.75, font = 2)
          }
        }
      }
    }

    graphics::box(lwd = 1.5)
  }

  # Panel 1: X -> Y (Row causes Column)
  graphics::par(mar = c(4, 4.5, 2, 0.5))
  draw_panel(mat_xy, mat_pval, mat_sig,
             expression(bold("Row") %->% bold("Column")),
             show_y_labels = TRUE)

  # Panel 2: Y -> X (Column causes Row) - use transposed matrix
  # Smaller left margin since no y-axis labels
  graphics::par(mar = c(4, 1, 2, 0.5))
  draw_panel(mat_yx, mat_pval_yx, mat_sig_yx,
             expression(bold("Column") %->% bold("Row")),
             show_y_labels = FALSE)

  # Color legend (only when gradient is enabled)
  if (gradient) {
    graphics::par(mar = c(5.5, 0.5, 3, 3))

    if (type == "significance") {
      graphics::image(1, c(0.25, 0.75), matrix(c(0, 1), ncol = 2),
                      col = colors, axes = FALSE, xlab = "", ylab = "")
      graphics::box()
      graphics::axis(4, at = c(0.25, 0.75), labels = c("Not\nSig.", "Sig."),
                     las = 1, cex.axis = 0.75, tick = FALSE, line = -0.5)
    } else {
      legend_seq <- seq(0, max_val, length.out = n_colors)
      graphics::image(1, legend_seq, matrix(legend_seq, ncol = length(legend_seq)),
                      col = colors, axes = FALSE, xlab = "", ylab = "")
      graphics::box()

      pretty_ticks <- pretty(c(0, max_val), n = 5)
      pretty_ticks <- pretty_ticks[pretty_ticks >= 0 & pretty_ticks <= max_val]
      graphics::axis(4, at = pretty_ticks, las = 1, cex.axis = 0.75)

      if (type == "pvalue") {
        graphics::mtext(expression(-log[10](p)), side = 4, line = 2, cex = 0.8)
        sig_threshold <- -log10(alpha)
        if (sig_threshold <= max_val) {
          graphics::abline(h = sig_threshold, col = "red", lty = 2, lwd = 1.5)
          graphics::text(1.5, sig_threshold, sprintf("p=%.2f", alpha),
                         col = "red", cex = 0.6, xpd = TRUE, pos = 4)
        }
      } else {
        graphics::mtext("F-stat", side = 4, line = 2, cex = 0.8)
      }
    }
  }

  invisible(x)
}
