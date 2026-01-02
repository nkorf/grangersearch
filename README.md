# grangersearch

An R package for Granger causality testing with tidyverse compatibility.

## Overview

`grangersearch` provides a simple, user-friendly interface for performing Granger causality tests on time series data. Key features include:

- **Pairwise testing**: Test whether one time series Granger-causes another
- **Exhaustive search**: Discover causal relationships among multiple variables
- **Lag selection**: Analyze sensitivity of results across different lag orders
- **Tidyverse integration**: Works seamlessly with pipes and data frames
- **Visualization**: Built-in plotting for causality matrices and lag analysis

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
set.seed(123)
x <- cumsum(rnorm(200))
y <- c(0, 0.7 * x[1:199]) + rnorm(200, sd = 0.5)

result <- granger_causality_test(x = x, y = y, lag = 2)
print(result)

# Tidyverse-style with data frames
library(tibble)
df <- tibble(x = x, y = y)
df |> granger_causality_test(x, y)

# Exhaustive search across multiple variables
data(Canada, package = "vars")
results <- granger_search(as.data.frame(Canada), lag = 2)
plot(results)  # Causality matrix visualization

# Lag selection analysis
lag_analysis <- granger_lag_select(df, x, y, lag = 1:8)
plot(lag_analysis)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `granger_causality_test()` | Test Granger causality between two time series |
| `granger_search()` | Exhaustive pairwise search across multiple variables |
| `granger_lag_select()` | Analyze results across different lag orders |
| `tidy()` / `glance()` | Broom-style tidying of results |

## Citation

If you use this package, please cite:

> Korfiatis, N. (2025). grangersearch: An R Package for Granger Causality Testing. arXiv preprint.

## Author

**Nikolaos Korfiatis**
Department of Informatics, Ionian University
Corfu, Greece
nkorf@ionio.gr

## License

MIT
