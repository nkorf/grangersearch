#' Perform a Granger causality test on two time series
#'
#' @param X a numeric vector or time series
#' @param Y a numeric vector or time series
#' @return a character string indicating the results of the test
#' @export

granger_causality_test <- function(X, Y) {
  # First, we will fit four VAR models:
  # model_X: only X as predictor
  # model_Y: only Y as predictor
  # model_XY: both X and Y as predictors
  # model_YX: both Y and X as predictors
  model_X <- VAR(X, p = 1)
  model_Y <- VAR(Y, p = 1)
  model_XY <- VAR(cbind(X, Y), p = 1)
  model_YX <- VAR(cbind(Y, X), p = 1)

  # Then we will perform the Granger causality test by comparing the AIC values of the four models
  granger_test_XY <- granger.test(model_XY, model_X, test = "F")
  granger_test_YX <- granger.test(model_YX, model_Y, test = "F")

  # The test results can be accessed through the 'p.value' attribute of the test object
  p_value_XY <- granger_test_XY$p.value
  p_value_YX <- granger_test_YX$p.value

  # We can conclude that X Granger causes Y if the p-value for model_XY is below a certain significance level
  if (p_value_XY < 0.05) {
    result <- "X Granger causes Y"
  } else {
    result <- "X does not Granger cause Y"
  }

  # We can conclude that Y Granger causes X if the p-value for model_YX is below a certain significance level
  if (p_value_YX < 0.05) {
    result <- paste(result, ", Y Granger causes X", sep = "")
  } else {
    result <- paste(result, ", Y does not Granger cause X", sep = "")
  }

  return(result)
}
