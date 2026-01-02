# grangersearch 0.1.0

Initial CRAN release.

## Features

* `granger_causality_test()` - Main function for performing Granger causality tests
  - Supports both vector input and tidyverse-style data frame input
  - Configurable lag order and significance level
  - Returns structured S3 object with p-values, test statistics, and conclusions

* Tidyverse integration
  - Pipe-compatible (`%>%` and `|>`)
  - Non-standard evaluation for column selection
  - `tidy()` method for broom-compatible tibble output
  - `glance()` method for model-level summary

* S3 methods
  - `print.granger_result()` for clean console output
  - `summary.granger_result()` for detailed results

* Example data
  - `example_causality` dataset with known causal relationship for testing and demonstration
