# grangersearch <img src="man/figures/logo.png" align="right" height="139" alt="grangersearch logo" />

An R package for exhaustive Granger causality testing with tidyverse integration.

## Overview

`grangersearch` provides a simple interface for performing Granger causality tests on time series data. The package wraps the `vars` infrastructure while providing a streamlined interface for exploratory causal analysis.

Key features include:

- **Exhaustive pairwise search**: Automatically discover Granger-causal relationships across multiple variables
- **Automatic lag optimization**: Systematic evaluation of multiple lag orders with visualization
- **Tidyverse compatibility**: Pipe operators (`|>`, `%>%`) and non-standard evaluation
- **Broom integration**: `tidy()` and `glance()` methods for structured output
- **Visualization**: Built-in plotting for causality matrices and lag selection analysis

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("nkorf/grangersearch")
```

## Quick Start

```r
library(grangersearch)

# Basic pairwise test
data(Canada, package = "vars")
result <- Canada |> granger_causality_test(e, U, lag = 2)
print(result)

# Get tidy results
tidy(result)

# Exhaustive search across multiple variables
search_results <- Canada |> granger_search(lag = 2)
plot(search_results)  # Causality matrix visualization

# Lag selection analysis
lag_analysis <- Canada |> granger_lag_select(e, U, lag = 1:8)
plot(lag_analysis)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `granger_causality_test()` | Test Granger causality between two time series |
| `granger_search()` | Exhaustive pairwise search across multiple variables |
| `granger_lag_select()` | Analyze results across different lag orders |
| `tidy()` / `glance()` | Broom-style tidying of results |

## Example Output

```
Granger Causality Test
======================

Observations: 84, Lag order: 2, Significance level: 0.050

e -> U: e Granger-causes U (p = 0.0000)
U -> e: U does not Granger-cause e (p = 0.2983)
```

## Citation

If you use this package, please cite:

> Korfiatis, N. (2025). grangersearch: An R Package for Exhaustive Granger Causality Testing with Tidyverse Integration. arXiv preprint. https://arxiv.org/abs/2601.01604

## Author

**Nikolaos Korfiatis**
Department of Informatics, Ionian University
Corfu, Greece
nkorf@ionio.gr

## License

MIT
