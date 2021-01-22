specs: Single-equation Penalized Error-Correction Selector
==========================================================

The R package `specs` implements the single-equation penalized
error-correction selector proposed in Smeekes and Wijler (2018a) as an
automated approach towards sparse single-equation cointegration
modelling. In addition, the package contains the dataset used in Smeekes
and Wijler (2018a) to predict Dutch unemployment rates based on Google
Trends.

Installation and Loading
------------------------

### Installation

The development version of the `specs` package can be installed from
GitHub using

    # install.package("devtools")
    devtools::install_github("wijler/specs")

### Load Package

After installation, the package can be loaded in the standard way:

    library(specs)

Data
----

Loading `specs` provides access to the dataset `Unempl_GT` which
contains the monthly unemployment rates (x1000) in the Netherlands and a
set of 87 related Google Trends aggregated to a monthly frequency. The
data covers the period from 1 January 2014 to December 2017 and is used
in the empirical application in Smeekes and Wijler (2018a).

Sparse single-equation cointegration modelling
----------------------------------------------

The package `specs` enables sparse and automated estimation of
single-equation cointegration models. To estimate a conditional
error-correction model (CECM) on the untransformed levels of a
collection of time series, use the function `specs`. Conversion of the
data matrix in levels to a CECM representation is performed
automatically within the function. For example,

    #library(specs)
    unemployment <- Unempl_GT[,1] #Extract the Dutch unemployment levels (x1,000)
    GT <- Unempl_GT[,2:11] #Select the first ten Google Trends
    my_specs <- specs(unemployment,GT,p=1) #Estimate specs
    my_coefs <- my_specs$gammas #store the coefficients

    y_d <- my_specs$y_d #Transformed dependent variable
    z_l <- my_specs$v #Transformed independent variables

estimates a regularized CECM on Dutch unemployment rates and ten Google
Trends by appropriately differencing and/or lagging the levels of the
time series, see Smeekes and Wijler (2018a, eq. 7) for the full model
specification. `my_coefs` contains 1,000 column, each column being the
solution for a unique combination of the group and individual penalty
(see the **penalty** section for more details).

The short-run dynamics in a CECM are modeled via the inclusion of lagged
differences of the data. By default, `specs` includes once lagged
differences, but the user is free to specify the desired number via the
input `p`. Continuing the above example,

    my_specs <- specs(unemployment,GT,p=0) #Estimate specs without any lagged differences

The inclusion of deterministic terms is often desired for correct model
specification. When it is believed that a constant and/or trend should
be included, it is advisable to estimate the model without imposing
regularization on these deterministic components. Accordingly, the input
`deterministics = c("constant","trend","both","none")` allows the user
to choose the deterministic specification without penalizing the
deterministic terms. For example,

    my_specs <- specs(unemployment,GT,deterministics = "both") #Estimate specs with a constant and trend included
    my_coefs <- my_specs$gammas #Store the new coefficients
    my_deterministics <- my_specs$thetas #Store the coefficients of the deterministic component

There are cases in which it may not be of interest to explicitly model
cointegration. In this instance, the lagged levels can be omitted from
the model by setting `ADL = TRUE`. This estimates a penalized
autoregressive distributed lag (ADL) model on the differenced data:

    my_adl <- specs(unemployment,GT,ADL=TRUE) #Estimate and ADL model
    my_coefs <- my_specs$gammas #Store the coefficients (smaller matrix than before)

The matrix of coefficients stored in `my_coefs` does not contain the
contribution of `z_l` anymore and, consequently, has a smaller row
dimension than before. `my_coefs` now contains only 100 columns, as
there is no more group penalty included. The reduced number of penalties
to estimate a solution far and some algorithmic change (see **Algorithm
and Implementation**) speed up the estimation procedure considerably.

Alternatively, the user may choose to pre-transform the data into the
form of a CECM/ADL, for example to save on computation time in rolling
window forecast exercises. In this case, the function `specs_tr` can be
used instead. This function operates entirely analogous to `specs`, with
the exception of requiring a differenced dependent variable (`y_d`), the
lagged levels of the time series (`z_l`) and the required differences of
the data (`w`) as inputs. Since `w` is directly provided to the
function, the option to set the lag length `p` is omitted. When
`ADL=TRUE`, the user may omit `z_l` as an input. For example:

    z <- cbind(unemployment,GT)
    y_d <- diff(unemployment) #Difference the dependent variable
    z_l <- GT[-nrow(GT),] #Lagged levels of the data
    w <- diff(GT) #Contemporaneous differences (corresponding to p=0)

    my_specs <- specs_tr(y_d,z_l,w) #Estimate a CECM on pre-transformed data
    my_adl <- specs_tr(y_d,NULL,w,ADL=TRUE) #Estimate an ADL model on pre-transformed data

Finally, a word on the naming of objects. Within this package the
coefficients corresponding to `z_l` are referred to as `delta`, those
corresponding to `w` as `pi`, and the numeric object that stacks both
`delta` and `pi` is referred to as `gamma`. Therefore, the solutions of
`specs` are referred to as `gammas`. The deterministic terms are passed
to the function outputs as `D`, with there coefficients being referred
to as `thetas`. As seen in the above examples, when `ADL=TRUE`, `delta`
is omitted from the numeric object `gamma` in the output. The naming of
objects here is congruent with Smeekes and Wijler (2018a), which may
serve as a helpful guideline for implementation.

**Penalty**

The functions `specs`, or equivalently `specs_tr`, estimate the model by
variants of penalized regression, customized to the error-correction
framework. In its most general form, specs penalizes each individual
coefficient via `lambda_i`, and adds a group penalty on `delta` via
`lambda_g`. Unless sequences of positive numbers for `lambda_i` and/or
`lambda_g` are supplied, sequences are generated automatically within
the function. A grid of 100 values is generated for `lambda_i` and a
grid of 10 values for `lambda_g`. The largest value in each grid
corresponds to the smallest value that sets all the coefficients that it
penalizes equal to zero (with the other penalty set equal to zero). The
smallest penalty in the grid is chosen as 1e-4 times the largest value
in that grid. As an important special case, the user may set
`lambda_g = 0`, in which case the function will estimate the model by
(weighted) *L*<sub>1</sub>-penalized regression, i.e. the (adaptive)
lasso:

    my_specs <- specs(unemployment,GT,lambda_g=0) #Estimate a CECM without group penalty

In practice, one typically requires a single choice of penalties that
provides the optimal solution for the model building exercise at hand.
To facilitate selection of such an optimal penalty, the functions
`specs_opt` and `specs_tr_opt`, being the equivalents of `specs` and
`specs_tr`, respectively, come with the added functionality of automated
selection of optimal values for `lambda_i` and `lambda_g`. The selection
criteria can be set via the input `rule`, with the possible choices
being `BIC`, `AIC` or time series cross-validation (`TSCV`). The first
two are information criteria in which the degrees of freedom is
approximated by the number of non-zero coefficients for a particular
solution, whereas the latter is a form of cross-validation that respects
the time series structure of the data. The implementation details for
TSCV can be found in Smeekes and Wijler (2018b, p. 411). The full matrix
of solutions, as well as the solution corresponding to the optimal
penalty choice are included in the output.

    my_specs <- specs_opt(unemployment,GT,rule="BIC") #Estimate a CECM with the optimal penalty chosen by BIC
    coefs_opt <- my_specs$gamma_opt #Extract the optimal coefficients

    my_adl <- specs_opt(unemployment,GT,rule="AIC",ADL=TRUE) #Estimate an ADL model with the optimal penalty chosen by BIC

    my_specs <- specs_opt(unemployment,GT,rule="TSCV",CV_cutoff=4/5) #Estimate a CECM with the optimal penalty chosen by TSCV
                                                                     #Training sample is 4/5 of the total sample

    my_specs <- specs_tr_opt(y_d,z_l,w,rule="BIC") #Estimate a CECM based on pre-transformed data, penalty chosen by BIC

**Weights**

Finally, specs can be estimated with the use of adaptively weighted
penalization. Automatically generated weighting schemes are available
via the input `weights = c("ridge","ols","none")`. The default option,
`"ridge"` constructs the weights via the use of initial estimates
obtained by ridge regression. In detail, the weights for
*δ*<sub>*i*</sub> and *π*<sub>*j*</sub> are constructed as
|*δ̂*<sub>*i*</sub>|<sup> − *k*<sub>*δ*</sub></sup> and
|*π̂*<sub>*j*</sub>|<sup> − *k*<sub>*π*</sub></sup>, respectively. The
penalty parameter for the ridge regression is automatically chosen by
TSCV. Alternative options are to automatically generate weights via
initial ols estimates (`weights = "ols"`) or to refrain from adaptive
weighting altogether (`weights = "none"`). Alternatively, it is also
possible to supply a sequence of positive weights directly. Finally, the
values for *k*<sub>*δ*</sub> and *k*<sub>*π*</sub> can be chosen by the
user via the equivalently named input options `k_delta` and `k_pi`. The
optimal values for these parameters are case-dependent, although some
theoretical guidance is provided in table 1 of Smeekes and Wijler
(2018a).

    my_specs <- specs(unemployment,GT,weights="ols",k_delta=2,k_pi=1) #Estimate specs with OLS and variable weight exponents

**Algorithm and Implementation**

The package `specs` combines accelerated generalized gradient descent
for the estimation of *δ* with coordinatewise descent for the estimation
of *π*. Since *δ* is penalized by both an *L*<sub>1</sub>- and
*L*<sub>2</sub>-penalty, its estimation fits into the framework of the
so-called sparse group lasso, for which numerous computational
procedures have been proposed. This package adopts the algorithm of
Simon et al. (2013), as the use of accelerated gradient descent via
Nesterov updates greatly improves computational time. However, since *π*
is regularized via an *L*<sub>1</sub>-penalty, and is separable from the
penalty on *δ*, the optimal solution for *π* is calculated via the
coordinate-wise descent procedure proposed in Friedman et al. (2009).
Essentially, `specs` iterates between optimizing for *δ* and *π*, where
within each iteration one of the aforementioned two algorithms is
repeated until numerical convergence. All calculations are performed in
C++, with the help of the Rcpp and Armadillo packages.

References
----------

-   Friedman, J., Hastie, T., and Tibshirani, R. (2009). glmnet: Lasso
    and elastic-net regularized generalized linear models. R package
    version, 1(4).
-   Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2013). A
    sparse-group lasso. Journal of computational and graphical
    statistics, 22(2), 231-245.
-   Smeekes, S., and Wijler, E. (2018a). An Automated Approach Towards
    Sparse Single-Equation Cointegration Modelling. arXiv preprint
    arXiv:1809.08889.
-   Smeekes, S., & Wijler, E. (2018b). Macroeconomic forecasting using
    penalized regression methods. International journal of forecasting,
    34(3), 408-430.
