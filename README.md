# grangeRsearch - R package

An R package for performing exhaustive search for granger causality for transactions time series. 
Given a `X` and `Y` timeseries as input, it performs a granger test to assertain causality. It considers the 
following manifestations: 

 * X: only X as predictor
 * Y: only Y as predictor
 * XY: both X and Y as predictors
 * YX: both Y and X as predictors

