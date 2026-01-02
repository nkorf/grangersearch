#' Example Time Series Data with Known Causal Relationship
#'
#' A dataset containing two time series where `cause_x` Granger-causes `effect_y`.
#' This data is useful for demonstrating and testing the Granger causality test.
#'
#' @format A data frame with 200 rows and 3 variables:
#' \describe{
#'   \item{time}{Integer. Time index from 1 to 200.}
#'   \item{cause_x}{Numeric. The "cause" time series, a random walk.}
#'   \item{effect_y}{Numeric. The "effect" time series, which depends on lagged values of cause_x.}
#' }
#'
#' @details
#' The data was generated with the following process:
#' \itemize{
#'   \item `cause_x` is a random walk: \eqn{X_t = X_{t-1} + \epsilon_t}
#'   \item `effect_y` depends on lagged `cause_x`: \eqn{Y_t = 0.7 \cdot X_{t-1} + \nu_t}
#' }
#' where \eqn{\epsilon_t \sim N(0, 1)} and \eqn{\nu_t \sim N(0, 0.5^2)}.
#'
#' When tested, `cause_x` should Granger-cause `effect_y`, but not vice versa.
#'
#' @examples
#' data(example_causality)
#'
#' # Test for Granger causality
#' result <- granger_causality_test(
#'   example_causality,
#'   cause_x,
#'   effect_y
#' )
#' print(result)
#'
#' @source Simulated data generated with seed 42.
"example_causality"
