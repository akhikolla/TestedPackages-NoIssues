# TSrepr 1.1.0
  * fixed repr_sma + added stopping rules with order parameter + unit tests on that
  * added new normalisation functions: arctan + Box-Cox + Yeo-Johnson
  * added repr_list function for handling multiple time series with different lengths
  
# TSrepr 1.0.4 2020/03/25

  * Fixed 0/0 case in forecasting accuracy measures
  * Refactor of data.table::melt cases in vignettes, because of data.table package changes
  * Added norm_*_params functions

# TSrepr 1.0.3 2019/05/31

  * New accuracy measure MSE (mean squared error) was added
  * Fixed SAX breaks (bins) length
  * Fixed some bad alignments in documentation
  * Added stopping criteria for forecasting accuracy measures, when real values and forecasts have a different lengths + tests for that
  * MAAPE is now consistent with MAPE and sMAPE, so it returns error in %

# TSrepr 1.0.2 2018/11/21

  * New accuracy measure MAAPE (mean arctangent absolute percentage error) was added
  * Added new references to vignettes
  * Added new references to documentation
  * Fixed some bad alignments in documentation

# TSrepr 1.0.1 2018/05/31

  * Created unit tests (by `testthat`) for all functions
  * Fixed vignette titles
  * New citation of package (CITATION file in \inst)
  * ORCID in DESCRIPTION

# TSrepr 1.0.0 2018/01/26

  * First CRAN release
