# borrowr

R package for estimating the population average treatment effect using a primary data source with borrowing from supplemental data sources.

To install from source and build vignettes:

devtools::install_github("jeffrey-boatman/borrowr", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), force = TRUE)

To do:

  - update documentation (changes included: adding argument for prior probability of exchangeability, updated gamma prior for bayesian linear model, ...)

