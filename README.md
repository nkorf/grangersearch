# grangeRsearch - R package

_Dr Nikolaos Korfiatis (n.korfiatis@uea.ac.uk)_
_Associate Professor in Business Analytics, University of East Anglia_ 
_Norwich Business School_

> This is a work in progress. Use at your own risk.

An R package for performing exhaustive search for granger causality for transactions time series. 
Given a `X` and `Y` timeseries as input, it performs a granger test to assertain causality. It considers the 
following manifestations: 

 * X: only X as predictor
 * Y: only Y as predictor
 * XY: both X and Y as predictors
 * YX: both Y and X as predictors

## Installation 

You can install the package using the `devtools::install_github("nkorf/grangersearch")`

