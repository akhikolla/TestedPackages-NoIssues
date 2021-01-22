## ----sim-data-----------------------------------------------------------------
# simulation constants
set.seed(467392)
n_obs <- 200
n_covars <- 3

# make some training data
x <- replicate(n_covars, rnorm(n_obs))
y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n_obs, mean = 0, sd = 0.2)

# make some testing data
test_x <- replicate(n_covars, rnorm(n_obs))
test_y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n_obs, mean = 0, sd = 0.2)

## ----sim-view-----------------------------------------------------------------
head(x)
head(y)

## -----------------------------------------------------------------------------
library(hal9001)

## ----fit-hal-glmnet-----------------------------------------------------------
hal_fit <- fit_hal(X = x, Y = y, fit_type = "glmnet")
hal_fit$times

## ----results-hal-glmnet-------------------------------------------------------
hal_fit

## ----fit-hal-reduced----------------------------------------------------------
hal_fit_reduced <- fit_hal(X = x, Y = y, fit_type = "lassi",
                           reduce_basis = 1/sqrt(length(y)))
hal_fit_reduced$times

## ----results-hal-reduced------------------------------------------------------
hal_fit_reduced

## ----eval-mse-----------------------------------------------------------------
# training sample prediction for HAL vs HAL9000
mse <- function(preds, y) {
    mean((preds - y)^2)
}

preds_hal <- predict(object = hal_fit, new_data = x)
mse_hal <- mse(preds = preds_hal, y = y)
mse_hal

## ----eval-oob-----------------------------------------------------------------
oob_hal <- predict(object = hal_fit, new_data = test_x)
oob_hal_mse <- mse(preds = oob_hal, y = test_y)
oob_hal_mse

