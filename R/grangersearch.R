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
#' @param type Character. Type of Granger causality to compute:
#'   \itemize{
#'     \item `"classic"` (default): Uses all lagged values from 1 to `lag` in the
#'       VAR model. This is the traditional Granger causality approach.
#'     \item `"constrained"`: Uses only the single value at the specified `lag`,
#'       not all values up to it. This approach has constant model complexity
#'       regardless of lag and tends to overfit less (Dimitrakopoulos, 2024).
#'   }
#' @param difference Logical. If TRUE, apply first-order differencing to both
#'   time series before analysis. This helps ensure stationarity, which is an
#'   assumption of Granger causality tests. Default is FALSE.
#'
#' @return An object of class `granger_result` containing:
#' \describe{
#'   \item{x_causes_y}{Logical. TRUE if X Granger-causes Y at the specified alpha level.}
#'   \item{y_causes_x}{Logical. TRUE if Y Granger-causes X at the specified alpha level.}
#'   \item{p_value_xy}{Numeric. P-value for the test of X causing Y.}
#'   \item{p_value_yx}{Numeric. P-value for the test of Y causing X.}
#'   \item{test_statistic_xy}{Numeric. Test statistic for X causing Y.}
#'   \item{test_statistic_yx}{Numeric. Test statistic for Y causing X.}
#'   \item{gc_strength_xy}{Numeric. Granger causality strength for X causing Y,
#'     computed as log(Var(univariate residuals) / Var(bivariate residuals)).
#'     Higher values indicate stronger predictive relationship.}
#'   \item{gc_strength_yx}{Numeric. Granger causality strength for Y causing X.}
#'   \item{lag}{Integer. The lag order used.}
#'   \item{alpha}{Numeric. The significance level used.}
#'   \item{test}{Character. The test type used.}
#'   \item{type}{Character. The type of Granger causality ("classic" or "constrained").}
#'   \item{difference}{Logical. Whether differencing was applied.}
#'   \item{n}{Integer. Number of observations (after differencing if applied).}
#'   \item{x_name}{Character. Name of the X variable.}
#'   \item{y_name}{Character. Name of the Y variable.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The Granger causality test is based on the idea that if X causes Y, then past
#' values of X should contain information that helps predict Y above and beyond
#' the information contained in past values of Y alone (Granger, 1969).
#'
#' For `type = "classic"`, this function fits Vector Autoregressive (VAR) models
#' using the \pkg{vars} package and performs F-tests to compare restricted and
#' unrestricted models. The test is performed in both directions to detect
#' unidirectional or bidirectional causality.
#'
#' For `type = "constrained"`, the function uses a simplified approach that only
#' considers the single lagged value at the specified lag (not all values from 1
#' to lag). This constrained approach has constant model complexity regardless of
#' the lag order and has been shown to overfit less than classic Granger causality,
#' especially for larger lag values. The overfitting behavior of classic GC models
#' is discussed in Shojaie & Fox (2022); the constrained formulation follows
#' Dimitrakopoulos (2024).
#'
#' The `gc_strength` values provide a continuous measure of Granger causality
#' magnitude, computed as log(Var(univariate) / Var(bivariate)). This formulation
#' follows Barrett et al. (2010). Higher values indicate that adding the predictor
#' variable substantially reduces prediction error variance.
#'
#' Note that Granger causality is a statistical concept based on prediction and
#' temporal precedence. It does not necessarily imply true causal mechanisms
#' (Granger, 1980).
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
#' # Access GC strength (continuous measure)
#' result$gc_strength_xy
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
#' # Use constrained Granger causality (less prone to overfitting)
#' df |> granger_causality_test(price, volume, type = "constrained", lag = 3)
#'
#' # Apply differencing for stationarity
#' df |> granger_causality_test(price, volume, difference = TRUE)
#'
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models
#' and cross-spectral methods. \emph{Econometrica}, 37(3), 424-438.
#'
#' Granger, C. W. J. (1980). Testing for causality: A personal viewpoint.
#' \emph{Journal of Economic Dynamics and Control}, 2, 329-352.
#'
#' Barrett, A. B., Barnett, L., & Seth, A. K. (2010). Multivariate Granger causality
#' and generalized variance. \emph{Physical Review E}, 81, 041907.
#'
#' Seth, A. K., Barrett, A. B., & Barnett, L. (2015). Granger causality analysis
#' in neuroscience and neuroimaging. \emph{Journal of Neuroscience}, 35(8), 3293-3297.
#'
#' Shojaie, A., & Fox, E. B. (2022). Granger causality: A review and recent advances.
#' \emph{Annual Review of Statistics and Its Application}, 9(1), 289-319.
#'
#' Dimitrakopoulos, P. S. (2024). Detecting Granger Causality. Master's Thesis,
#' Eindhoven University of Technology.
#'
#' @seealso
#' \code{\link[vars]{VAR}} for the underlying VAR model,
#' \code{\link[vars]{causality}} for an alternative implementation,
#' [tidy.granger_result()] for tidying results.
#'
#' @importFrom rlang enquo eval_tidy as_label
#' @importFrom stats complete.cases lm var pf
#' @export
granger_causality_test <- function(.data = NULL, x, y, lag = 1, alpha = 0.05,
                                    test = "F", type = c("classic", "constrained"),
                                    difference = FALSE) {
  type <- match.arg(type)
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

  if (!is.logical(difference) || length(difference) != 1) {
    stop("`difference` must be TRUE or FALSE.", call. = FALSE)
  }

  if (any(is.na(x_vec)) || any(is.nan(x_vec)) || any(is.infinite(x_vec))) {
    stop(sprintf("`%s` contains NA, NaN, or infinite values.", x_name), call. = FALSE)
  }
  if (any(is.na(y_vec)) || any(is.nan(y_vec)) || any(is.infinite(y_vec))) {
    stop(sprintf("`%s` contains NA, NaN, or infinite values.", y_name), call. = FALSE)
  }

  # Apply differencing if requested (for stationarity)
  if (difference) {
    x_vec <- diff(x_vec)
    y_vec <- diff(y_vec)
  }

  n <- length(x_vec)
  min_length <- 2 * lag + 2
  if (n < min_length) {
    stop(sprintf("Time series too short. Need at least %d observations for lag = %d%s.",
                 min_length, lag,
                 if (difference) " (after differencing)" else ""), call. = FALSE)
  }

  # Initialize results
  p_value_xy <- NA_real_
  p_value_yx <- NA_real_
  test_stat_xy <- NA_real_
  test_stat_yx <- NA_real_
  gc_strength_xy <- NA_real_
  gc_strength_yx <- NA_real_

  if (type == "classic") {
    # Classic Granger causality using VAR models
    data_matrix <- cbind(X = x_vec, Y = y_vec)

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

      # Calculate GC strength from residual variances
      # For X -> Y: compare Y equation with and without X lags
      resid_y <- stats::residuals(var_model)[, "Y"]
      var_bivariate_y <- stats::var(resid_y)

      # Fit univariate AR model for Y
      y_lagged <- stats::embed(y_vec, lag + 1)
      y_response <- y_lagged[, 1]
      y_predictors <- y_lagged[, -1, drop = FALSE]
      uni_model_y <- stats::lm(y_response ~ y_predictors)
      var_univariate_y <- stats::var(stats::residuals(uni_model_y))
      gc_strength_xy <- log(var_univariate_y / var_bivariate_y)

      # For Y -> X: compare X equation with and without Y lags
      resid_x <- stats::residuals(var_model)[, "X"]
      var_bivariate_x <- stats::var(resid_x)

      x_lagged <- stats::embed(x_vec, lag + 1)
      x_response <- x_lagged[, 1]
      x_predictors <- x_lagged[, -1, drop = FALSE]
      uni_model_x <- stats::lm(x_response ~ x_predictors)
      var_univariate_x <- stats::var(stats::residuals(uni_model_x))
      gc_strength_yx <- log(var_univariate_x / var_bivariate_x)

    }, error = function(e) {
      stop(sprintf("VAR model fitting failed: %s", e$message), call. = FALSE)
    })

  } else {
    # Constrained Granger causality: only use value at lag q (not all 1:q)
    # This has constant model complexity regardless of lag

    # Create lagged variables at exactly lag q
    effective_n <- n - lag
    y_response <- y_vec[(lag + 1):n]
    y_lag_q <- y_vec[1:effective_n]
    x_lag_q <- x_vec[1:effective_n]

    # Test X -> Y (does X help predict Y?)
    # Univariate model: Y_t ~ Y_{t-q}
    uni_model_xy <- stats::lm(y_response ~ y_lag_q)
    resid_uni_xy <- stats::residuals(uni_model_xy)
    var_uni_xy <- stats::var(resid_uni_xy)

    # Bivariate model: Y_t ~ Y_{t-q} + X_{t-q}
    bi_model_xy <- stats::lm(y_response ~ y_lag_q + x_lag_q)
    resid_bi_xy <- stats::residuals(bi_model_xy)
    var_bi_xy <- stats::var(resid_bi_xy)

    # GC strength: log ratio of variances
    gc_strength_xy <- log(var_uni_xy / var_bi_xy)

    # F-test for X -> Y
    # Compare nested models: does adding x_lag_q improve the model?
    anova_xy <- stats::anova(uni_model_xy, bi_model_xy)
    test_stat_xy <- anova_xy$F[2]
    p_value_xy <- anova_xy$`Pr(>F)`[2]

    # Test Y -> X (does Y help predict X?)
    x_response <- x_vec[(lag + 1):n]

    # Univariate model: X_t ~ X_{t-q}
    uni_model_yx <- stats::lm(x_response ~ x_lag_q)
    resid_uni_yx <- stats::residuals(uni_model_yx)
    var_uni_yx <- stats::var(resid_uni_yx)

    # Bivariate model: X_t ~ X_{t-q} + Y_{t-q}
    bi_model_yx <- stats::lm(x_response ~ x_lag_q + y_lag_q)
    resid_bi_yx <- stats::residuals(bi_model_yx)
    var_bi_yx <- stats::var(resid_bi_yx)

    # GC strength for Y -> X
    gc_strength_yx <- log(var_uni_yx / var_bi_yx)

    # F-test for Y -> X
    anova_yx <- stats::anova(uni_model_yx, bi_model_yx)
    test_stat_yx <- anova_yx$F[2]
    p_value_yx <- anova_yx$`Pr(>F)`[2]
  }

  # Build result object
  result <- list(
    x_causes_y = p_value_xy < alpha,
    y_causes_x = p_value_yx < alpha,
    p_value_xy = as.numeric(p_value_xy),
    p_value_yx = as.numeric(p_value_yx),
    test_statistic_xy = as.numeric(test_stat_xy),
    test_statistic_yx = as.numeric(test_stat_yx),
    gc_strength_xy = as.numeric(gc_strength_xy),
    gc_strength_yx = as.numeric(gc_strength_yx),
    lag = lag,
    alpha = alpha,
    test = test,
    type = type,
    difference = difference,
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

  type_label <- if (!is.null(x$type) && x$type == "constrained") "constrained" else "classic"
  diff_label <- if (!is.null(x$difference) && x$difference) ", differenced" else ""

  cat(sprintf("Type: %s, Observations: %d, Lag: %d, Alpha: %.3f%s\n\n",
              type_label, x$n, x$lag, x$alpha, diff_label))

  # X -> Y result
  xy_result <- if (x$x_causes_y) {
    sprintf("%s Granger-causes %s", x$x_name, x$y_name)
  } else {
    sprintf("%s does not Granger-cause %s", x$x_name, x$y_name)
  }
  gc_xy <- if (!is.null(x$gc_strength_xy) && !is.na(x$gc_strength_xy)) {
    sprintf(", GC = %.4f", x$gc_strength_xy)
  } else ""
  cat(sprintf("%s -> %s: %s (p = %.4f%s)\n", x$x_name, x$y_name, xy_result, x$p_value_xy, gc_xy))

  # Y -> X result
  yx_result <- if (x$y_causes_x) {
    sprintf("%s Granger-causes %s", x$y_name, x$x_name)
  } else {
    sprintf("%s does not Granger-cause %s", x$y_name, x$x_name)
  }
  gc_yx <- if (!is.null(x$gc_strength_yx) && !is.na(x$gc_strength_yx)) {
    sprintf(", GC = %.4f", x$gc_strength_yx)
  } else ""
  cat(sprintf("%s -> %s: %s (p = %.4f%s)\n", x$y_name, x$x_name, yx_result, x$p_value_yx, gc_yx))

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

  type_label <- if (!is.null(object$type) && object$type == "constrained") "constrained" else "classic"

  cat(sprintf("Variables: %s, %s\n", object$x_name, object$y_name))
  cat(sprintf("Number of observations: %d\n", object$n))
  cat(sprintf("Lag order: %d\n", object$lag))
  cat(sprintf("GC type: %s\n", type_label))
  cat(sprintf("Test type: %s-test\n", object$test))
  cat(sprintf("Significance level (alpha): %.3f\n", object$alpha))
  if (!is.null(object$difference) && object$difference) {
    cat("Differencing: applied (first-order)\n")
  }
  cat("\n")

  cat("Results:\n")
  cat("--------\n\n")

  cat(sprintf("Test: %s Granger-causes %s\n", object$x_name, object$y_name))
  cat(sprintf("  Test statistic: %.4f\n", object$test_statistic_xy))
  cat(sprintf("  P-value: %.6f\n", object$p_value_xy))
  if (!is.null(object$gc_strength_xy) && !is.na(object$gc_strength_xy)) {
    cat(sprintf("  GC strength: %.4f\n", object$gc_strength_xy))
  }
  cat(sprintf("  Conclusion: %s\n\n",
              if (object$x_causes_y) {
                sprintf("REJECT null (%s causes %s)", object$x_name, object$y_name)
              } else {
                sprintf("FAIL TO REJECT null (%s does not cause %s)", object$x_name, object$y_name)
              }))

  cat(sprintf("Test: %s Granger-causes %s\n", object$y_name, object$x_name))
  cat(sprintf("  Test statistic: %.4f\n", object$test_statistic_yx))
  cat(sprintf("  P-value: %.6f\n", object$p_value_yx))
  if (!is.null(object$gc_strength_yx) && !is.na(object$gc_strength_yx)) {
    cat(sprintf("  GC strength: %.4f\n", object$gc_strength_yx))
  }
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
#'   \item{gc_strength}{Numeric. The Granger causality strength (log variance ratio).}
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
  gc_xy <- if (!is.null(x$gc_strength_xy)) x$gc_strength_xy else NA_real_
  gc_yx <- if (!is.null(x$gc_strength_yx)) x$gc_strength_yx else NA_real_

  tibble::tibble(
    direction = c(
      sprintf("%s -> %s", x$x_name, x$y_name),
      sprintf("%s -> %s", x$y_name, x$x_name)
    ),
    cause = c(x$x_name, x$y_name),
    effect = c(x$y_name, x$x_name),
    statistic = c(x$test_statistic_xy, x$test_statistic_yx),
    p.value = c(x$p_value_xy, x$p_value_yx),
    gc_strength = c(gc_xy, gc_yx),
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
#'   \item{type}{Character. GC type used ("classic" or "constrained").}
#'   \item{difference}{Logical. Whether differencing was applied.}
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
  type_val <- if (!is.null(x$type)) x$type else "classic"
  diff_val <- if (!is.null(x$difference)) x$difference else FALSE

  tibble::tibble(
    nobs = x$n,
    lag = x$lag,
    alpha = x$alpha,
    test = x$test,
    type = type_val,
    difference = diff_val,
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
#' @param type Character. Type of Granger causality: "classic" (default) or
#'   "constrained". See [granger_causality_test()] for details.
#' @param difference Logical. If TRUE, apply first-order differencing before
#'   analysis. Default is FALSE.
#' @param include_insignificant Logical. If FALSE (default), only return
#'   significant causal relationships. If TRUE, return all pairwise results.
#'
#' @return A tibble with one row per directed pair tested, containing:
#' \describe{
#'   \item{cause}{Character. The potential cause variable name.}
#'   \item{effect}{Character. The potential effect variable name.}
#'   \item{statistic}{Numeric. The F-test statistic.}
#'   \item{p.value}{Numeric. The p-value of the test.}
#'   \item{gc_strength}{Numeric. The Granger causality strength (log variance ratio).}
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
#' # Use constrained Granger causality (less prone to overfitting)
#' granger_search(df, type = "constrained")
#'
#' # Apply differencing for stationarity
#' granger_search(df, difference = TRUE)
#'
#' @seealso [granger_causality_test()] for testing a single pair.
#'
#' @importFrom tibble tibble
#' @importFrom rlang enquos eval_tidy as_label
#' @export
granger_search <- function(.data, ..., lag = 1, alpha = 0.05, test = "F",
                           type = c("classic", "constrained"),
                           difference = FALSE,
                           include_insignificant = FALSE) {
  type <- match.arg(type)

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  # Validate lag parameter
  if (!is.numeric(lag) || any(lag < 1) || any(lag != as.integer(lag))) {
    stop("`lag` must be positive integer(s).", call. = FALSE)
  }
  lag <- as.integer(lag)

  if (!is.logical(difference) || length(difference) != 1) {
    stop("`difference` must be TRUE or FALSE.", call. = FALSE)
  }

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
          lag = p, alpha = alpha, test = test,
          type = type, difference = difference
        )

        if (result$p_value_xy < best_p) {
          best_p <- result$p_value_xy
          best_result <- tibble::tibble(
            cause = cause_var,
            effect = effect_var,
            statistic = result$test_statistic_xy,
            p.value = result$p_value_xy,
            gc_strength = result$gc_strength_xy,
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
        gc_strength = NA_real_,
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
  attr(output, "type") <- type
  attr(output, "difference") <- difference

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
#' @param type Character. Type of Granger causality: "classic" (default) or
#'   "constrained". See [granger_causality_test()] for details.
#' @param difference Logical. If TRUE, apply first-order differencing before
#'   analysis. Default is FALSE.
#'
#' @return A tibble with one row per (cause, effect, lag) combination:
#' \describe{
#'   \item{cause}{Character. The potential cause variable.}
#'   \item{effect}{Character. The potential effect variable.}
#'   \item{lag}{Integer. The lag order tested.}
#'   \item{statistic}{Numeric. The F-test statistic.}
#'   \item{p.value}{Numeric. The p-value.}
#'   \item{gc_strength}{Numeric. The Granger causality strength (log variance ratio).}
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
granger_lag_select <- function(.data, ..., lag = 1:4, alpha = 0.05, test = "F",
                                type = c("classic", "constrained"),
                                difference = FALSE) {
  type <- match.arg(type)

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(lag) || any(lag < 1) || any(lag != as.integer(lag))) {
    stop("`lag` must be positive integer(s).", call. = FALSE)
  }
  lag <- as.integer(lag)

  if (!is.logical(difference) || length(difference) != 1) {
    stop("`difference` must be TRUE or FALSE.", call. = FALSE)
  }

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
          lag = p, alpha = alpha, test = test,
          type = type, difference = difference
        )

        results[[length(results) + 1]] <- tibble::tibble(
          cause = cause_var,
          effect = effect_var,
          lag = p,
          statistic = result$test_statistic_xy,
          p.value = result$p_value_xy,
          gc_strength = result$gc_strength_xy,
          significant = result$x_causes_y
        )
      }, error = function(e) {
        results[[length(results) + 1]] <<- tibble::tibble(
          cause = cause_var,
          effect = effect_var,
          lag = p,
          statistic = NA_real_,
          p.value = NA_real_,
          gc_strength = NA_real_,
          significant = NA
        )
      })
    }
  }

  output <- do.call(rbind, results)
  attr(output, "alpha") <- alpha
  attr(output, "type") <- type
  attr(output, "difference") <- difference
  class(output) <- c("granger_lag_select", class(output))
  output
}


#' Plot Lag Selection Results
#'
#' Creates a visualization of p-values across different lag orders for
#' Granger causality tests.
#'
#' @param x A `granger_lag_select` object from [granger_lag_select()].
#' @param monochrome Logical. If TRUE, uses black lines with different line types
#'   and point symbols instead of colors. Suitable for publications. Default FALSE.
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
plot.granger_lag_select <- function(x, monochrome = FALSE, ...) {
  alpha <- attr(x, "alpha") %||% 0.05

  # Get unique directions
  x$direction <- paste(x$cause, "->", x$effect)
  directions <- unique(x$direction)
  n_dir <- length(directions)

  # Set up colors and line types
  if (monochrome) {
    colors <- rep("black", n_dir)
    ltys <- rep(1:6, length.out = n_dir)
    pchs <- c(1, 2, 0, 5, 6, 4)[rep(1:6, length.out = n_dir)]
  } else {
    colors <- grDevices::rainbow(n_dir)
    ltys <- rep(1, n_dir)
    pchs <- rep(19, n_dir)
  }

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
    graphics::lines(subset_data$lag, subset_data$p.value, col = colors[i], lwd = 2, lty = ltys[i])
    graphics::points(subset_data$lag, subset_data$p.value, col = colors[i], pch = pchs[i])
  }

  # Add legend
  graphics::legend(
    "topright",
    legend = directions,
    col = colors,
    lwd = 2,
    lty = ltys,
    pch = pchs,
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
#' @param monochrome Logical. If TRUE, uses grayscale palette instead of blue.
#'   Suitable for publications. Default FALSE.
#' @param show_metric Character. Additional metric to show in parentheses:
#'   "none" (default), "pvalue", or "strength" (GC strength).
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
                                        show_values = TRUE, gradient = TRUE,
                                        monochrome = FALSE,
                                        show_metric = c("none", "pvalue", "strength"), ...) {
  type <- match.arg(type)
  show_metric <- match.arg(show_metric)
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
  mat_strength <- matrix(NA, nrow = n_vars, ncol = n_vars,
                         dimnames = list(all_vars, all_vars))

  # Fill matrices
  for (i in seq_len(nrow(x))) {
    cause <- x$cause[i]
    effect <- x$effect[i]

    mat_pval[cause, effect] <- x$p.value[i]
    mat_sig[cause, effect] <- as.numeric(x$significant[i])
    if ("gc_strength" %in% names(x)) {
      mat_strength[cause, effect] <- x$gc_strength[i]
    }

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
  mat_strength_yx <- t(mat_strength)

  # Calculate max value for color scaling (same scale for both)
  max_val <- max(mat_xy, na.rm = TRUE)
  if (is.na(max_val) || max_val <= 0) max_val <- 1

  # Set up layout: two panels (+ optional color legend at bottom)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (gradient) {
    # Layout: horizontal legend on top, two heatmaps on bottom row
    graphics::layout(matrix(c(3, 3, 1, 2), nrow = 2, byrow = TRUE),
                     widths = c(1.3, 1), heights = c(0.5, 10))
  } else {
    # No legend needed when gradient is off
    # Left panel wider to accommodate y-axis labels
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.3, 1))
  }

  # Color palette
  if (!gradient || type == "significance") {
    if (monochrome) {
      colors <- c("gray90", "gray30")
    } else {
      colors <- c("gray90", "steelblue3")
    }
    breaks <- c(-0.5, 0.5, 1.5)
    n_colors <- 2
  } else {
    n_colors <- 100
    if (monochrome) {
      colors <- grDevices::colorRampPalette(c("gray95", "gray70", "gray40", "gray10"))(n_colors)
    } else {
      colors <- grDevices::colorRampPalette(c("gray95", "lightblue", "steelblue", "darkblue"))(n_colors)
    }
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
              label <- if (!is.na(sig_val) && sig_val > 0.5) "*" else ""
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

  # Color legend (only when gradient is enabled) - horizontal at top
  if (gradient) {
    graphics::par(mar = c(0.1, 6, 1.5, 6))

    if (type == "significance") {
      graphics::image(c(0.25, 0.75), 1, matrix(c(0, 1), nrow = 2),
                      col = colors, axes = FALSE, xlab = "", ylab = "")
      graphics::axis(3, at = c(0.25, 0.75), labels = c("Not Sig.", "Sig."),
                     cex.axis = 0.6, tick = FALSE, line = -0.5)
    } else {
      legend_seq <- seq(0, max_val, length.out = n_colors)
      graphics::image(legend_seq, 1, matrix(legend_seq, nrow = length(legend_seq)),
                      col = colors, axes = FALSE, xlab = "", ylab = "")

      pretty_ticks <- pretty(c(0, max_val), n = 5)
      pretty_ticks <- pretty_ticks[pretty_ticks >= 0 & pretty_ticks <= max_val]
      graphics::axis(3, at = pretty_ticks, cex.axis = 0.6)

      if (type == "pvalue") {
        graphics::mtext(expression(-log[10](p)), side = 3, line = 0.8, cex = 0.65)
        sig_threshold <- -log10(alpha)
        if (sig_threshold <= max_val) {
          graphics::abline(v = sig_threshold, col = "red", lty = 2, lwd = 1.5)
        }
      } else {
        graphics::mtext("F-stat", side = 3, line = 0.8, cex = 0.65)
      }
    }
  }

  invisible(x)
}


#' Distribution Analysis of Granger Causality
#'
#' Computes Granger causality for all pairwise combinations and returns
#' detailed distribution information. This is useful for understanding the
#' overall pattern of causal relationships in a dataset.
#'
#' @param .data A data frame or tibble containing the time series variables.
#' @param ... <[`tidy-select`][dplyr::dplyr_tidy_select]> Columns to include.
#'   If empty, all numeric columns are used.
#' @param lag Integer or integer vector. The lag order(s) to test. If a vector,
#'   results are returned for each lag separately. Default is 1.
#' @param type Character. Type of Granger causality: "classic" (default) or
#'   "constrained". See [granger_causality_test()] for details.
#' @param difference Logical. If TRUE, apply first-order differencing before
#'   analysis. Default is FALSE.
#'
#' @return An object of class `granger_distribution` containing:
#' \describe{
#'   \item{data}{A tibble with all pairwise GC results including gc_strength values.}
#'   \item{summary}{A tibble with summary statistics for each lag.}
#'   \item{lag}{The lag(s) used.}
#'   \item{type}{The type of GC computed.}
#'   \item{difference}{Whether differencing was applied.}
#'   \item{n_vars}{Number of variables analyzed.}
#'   \item{n_pairs}{Number of directed pairs tested.}
#' }
#'
#' @details
#' This function is designed for exploratory analysis of Granger causality
#' distributions across a dataset. It computes the GC strength (log variance
#' ratio) for all directed pairs and provides summary statistics.
#'
#' The distribution analysis helps understand:
#' \itemize{
#'   \item The overall spread of GC values in the dataset
#'   \item How GC distributions change across different lags
#'   \item Whether there are outliers with unusually high GC values
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' df <- data.frame(
#'   A = cumsum(rnorm(n)),
#'   B = cumsum(rnorm(n)),
#'   C = cumsum(rnorm(n))
#' )
#' # Add causal structure
#' df$B <- c(0, 0.7 * df$A[1:(n-1)]) + rnorm(n, sd = 0.5)
#'
#' # Analyze GC distribution
#' dist <- granger_distribution(df)
#' print(dist)
#'
#' # Visualize the distribution
#' plot(dist)
#'
#' # Compare classic vs constrained across lags
#' dist_classic <- granger_distribution(df, lag = 1:5, type = "classic")
#' dist_constrained <- granger_distribution(df, lag = 1:5, type = "constrained")
#'
#' @seealso [granger_search()] for finding significant relationships,
#'   [plot.granger_distribution()] for visualization.
#'
#' @importFrom tibble tibble
#' @importFrom stats quantile median
#' @export
granger_distribution <- function(.data, ..., lag = 1,
                                  type = c("classic", "constrained"),
                                  difference = FALSE) {
  type <- match.arg(type)

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(lag) || any(lag < 1) || any(lag != as.integer(lag))) {
    stop("`lag` must be positive integer(s).", call. = FALSE)
  }
  lag <- as.integer(lag)

  if (!is.logical(difference) || length(difference) != 1) {
    stop("`difference` must be TRUE or FALSE.", call. = FALSE)
  }

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

  n_vars <- length(var_names)

  # Generate all pairs
  pairs <- expand.grid(cause = var_names, effect = var_names,
                       stringsAsFactors = FALSE)
  pairs <- pairs[pairs$cause != pairs$effect, ]
  n_pairs <- nrow(pairs)

  # Compute GC for all pairs at all lags
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
          lag = p, alpha = 0.05, test = "F",
          type = type, difference = difference
        )

        results[[length(results) + 1]] <- tibble::tibble(
          cause = cause_var,
          effect = effect_var,
          lag = p,
          gc_strength = result$gc_strength_xy,
          p.value = result$p_value_xy,
          statistic = result$test_statistic_xy
        )
      }, error = function(e) {
        results[[length(results) + 1]] <<- tibble::tibble(
          cause = cause_var,
          effect = effect_var,
          lag = p,
          gc_strength = NA_real_,
          p.value = NA_real_,
          statistic = NA_real_
        )
      })
    }
  }

  data <- do.call(rbind, results)

  # Compute summary statistics by lag
  summary_list <- lapply(lag, function(p) {
    subset_data <- data[data$lag == p & !is.na(data$gc_strength), ]
    gc_vals <- subset_data$gc_strength

    if (length(gc_vals) > 0) {
      tibble::tibble(
        lag = p,
        n = length(gc_vals),
        mean = mean(gc_vals),
        median = stats::median(gc_vals),
        sd = stats::sd(gc_vals),
        min = min(gc_vals),
        max = max(gc_vals),
        q25 = stats::quantile(gc_vals, 0.25),
        q75 = stats::quantile(gc_vals, 0.75)
      )
    } else {
      tibble::tibble(
        lag = p,
        n = 0L,
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        min = NA_real_,
        max = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_
      )
    }
  })
  summary_data <- do.call(rbind, summary_list)

  result <- list(
    data = data,
    summary = summary_data,
    lag = lag,
    type = type,
    difference = difference,
    n_vars = n_vars,
    n_pairs = n_pairs
  )

  class(result) <- "granger_distribution"
  result
}


#' Print Method for granger_distribution Objects
#'
#' @param x A `granger_distribution` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.granger_distribution <- function(x, ...) {
  cat("\nGranger Causality Distribution Analysis\n")
  cat("========================================\n\n")

  cat(sprintf("Type: %s\n", x$type))
  cat(sprintf("Variables: %d, Directed pairs: %d\n", x$n_vars, x$n_pairs))
  if (x$difference) {
    cat("Differencing: applied\n")
  }
  cat(sprintf("Lag(s) tested: %s\n\n", paste(x$lag, collapse = ", ")))

  cat("Summary statistics (GC strength):\n")
  cat("---------------------------------\n")
  print(x$summary)

  cat("\n")
  invisible(x)
}


#' Plot Granger Causality Distribution
#'
#' Creates visualizations of the Granger causality strength distribution.
#'
#' @param x A `granger_distribution` object from [granger_distribution()].
#' @param type Character. Type of plot: "histogram" (default), "density", or "violin".
#' @param monochrome Logical. If TRUE, uses black/gray instead of colors.
#'   Suitable for publications. Default FALSE.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @details
#' For multiple lags:
#' \itemize{
#'   \item "histogram": Faceted histograms by lag
#'   \item "density": Overlaid density curves colored by lag
#'   \item "violin": Violin plots comparing distributions across lags
#' }
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   A = cumsum(rnorm(100)),
#'   B = cumsum(rnorm(100)),
#'   C = cumsum(rnorm(100))
#' )
#' df$B <- c(0, 0.7 * df$A[1:99]) + rnorm(100, sd = 0.5)
#'
#' dist <- granger_distribution(df, lag = 1:3)
#' plot(dist)
#' plot(dist, type = "density")
#' plot(dist, type = "violin")
#'
#' @export
plot.granger_distribution <- function(x, type = c("histogram", "density", "violin"),
                                       monochrome = FALSE, ...) {
  type <- match.arg(type)

  data <- x$data[!is.na(x$data$gc_strength), ]
  n_lags <- length(x$lag)

  if (nrow(data) == 0) {
    message("No valid GC values to plot.")
    return(invisible(x))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (type == "histogram") {
    if (n_lags > 1) {
      n_cols <- min(3, n_lags)
      n_rows <- ceiling(n_lags / n_cols)
      graphics::par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
    }

    hist_col <- if (monochrome) "gray60" else "steelblue"
    hist_border <- if (monochrome) "gray30" else "white"
    mean_col <- if (monochrome) "black" else "red"

    for (p in x$lag) {
      subset_data <- data[data$lag == p, ]
      gc_vals <- subset_data$gc_strength

      if (length(gc_vals) > 0) {
        graphics::hist(
          gc_vals,
          main = sprintf("Lag = %d (%s)", p, x$type),
          xlab = "GC Strength (log variance ratio)",
          col = hist_col,
          border = hist_border,
          breaks = "Sturges"
        )
        graphics::abline(v = mean(gc_vals), col = mean_col, lwd = 2, lty = 2)
        graphics::abline(v = 0, col = "gray50", lwd = 1, lty = 3)
      }
    }

  } else if (type == "density") {
    # Set up colors and line types
    if (monochrome) {
      colors <- rep("black", n_lags)
      ltys <- rep(1:6, length.out = n_lags)
    } else {
      colors <- grDevices::rainbow(n_lags)
      ltys <- rep(1, n_lags)
    }

    # Get overall range
    all_gc <- data$gc_strength
    x_range <- range(all_gc)

    # Calculate densities first (without plotting)
    max_density <- 0
    densities <- list()

    for (i in seq_along(x$lag)) {
      p <- x$lag[i]
      subset_data <- data[data$lag == p, ]
      gc_vals <- subset_data$gc_strength

      if (length(gc_vals) > 1) {
        d <- stats::density(gc_vals)
        densities[[i]] <- d
        max_density <- max(max_density, max(d$y))
      }
    }

    # Single plot with correct y-axis
    graphics::plot(
      x_range, c(0, max_density * 1.1),
      type = "n",
      xlab = "GC Strength (log variance ratio)",
      ylab = "Density",
      main = sprintf("GC Distribution by Lag (%s)", x$type)
    )

    for (i in seq_along(densities)) {
      if (!is.null(densities[[i]])) {
        graphics::lines(densities[[i]], col = colors[i], lwd = 2, lty = ltys[i])
      }
    }

    graphics::abline(v = 0, col = "gray50", lwd = 1, lty = 3)

    graphics::legend(
      "topright",
      legend = paste("Lag", x$lag),
      col = colors,
      lwd = 2,
      lty = ltys,
      cex = 0.8,
      bg = "white"
    )

  } else if (type == "violin") {
    # Simple boxplot-based visualization (violin requires additional packages)
    data$lag_factor <- factor(data$lag, levels = x$lag)

    box_col <- if (monochrome) "gray60" else "steelblue"
    box_border <- if (monochrome) "black" else "darkblue"

    graphics::boxplot(
      gc_strength ~ lag_factor,
      data = data,
      main = sprintf("GC Distribution by Lag (%s)", x$type),
      xlab = "Lag",
      ylab = "GC Strength (log variance ratio)",
      col = box_col,
      border = box_border
    )
    graphics::abline(h = 0, col = "gray50", lwd = 1, lty = 3)
  }

  invisible(x)
}
