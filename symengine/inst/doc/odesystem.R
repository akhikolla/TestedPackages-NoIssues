## -----------------------------------------------------------------------------
library(symengine)

## -----------------------------------------------------------------------------
x <- Symbol("x")
a <- 3
eq <- dxdt(x) == a/(x + 1)
eq

## -----------------------------------------------------------------------------
sigma <- 10
rho <- 28
beta <- 8/3
use_vars(x, y, z)

## -----------------------------------------------------------------------------
lorenz_sys <- list(
    dxdt(x) == sigma * (y - x),
    dxdt(y) == (rho - z) * x - y,
    dxdt(z) == - beta * z + x * y
)
lorenz_sys <- ODESystem(lorenz_sys, method = "rk5_i")

## -----------------------------------------------------------------------------
res <- predict(lorenz_sys, init = c(x=1, y=1, z=1),
               duration = 100, step_size = 0.001, start = 0)
head(res)

## ----fig.height=5, fig.width=5------------------------------------------------
plot(res[, c(2, 4)], type = 'l', col = "steelblue", main = "Lorenz Attractor")

## -----------------------------------------------------------------------------
use_vars(x, y)
vdp_sys <- ODESystem(
    dxdt(x) == y,
    dxdt(y) == 2 * (1 - x * x) * y - x,
    method = "bsd" # Bulirsch-Stoer
)
res <- predict(vdp_sys, init = rep(1e-4, 2), duration = 100, step_size = 0.01)

## ----fig.height=5, fig.width=5------------------------------------------------
oldpar <- par(mfrow = c(2, 2), mar = rep(0.5, 4), oma = rep(5, 4), xpd = NA)
make.plot <- function(xy, xlab = NA, ylab = NA)
  plot(xy, col = "steelblue", lwd = 2, type = "l",
       axes = FALSE, xlab = xlab, ylab = ylab)
plot.new()
make.plot(res[, c(3, 1)]); axis(3); axis(4)
make.plot(res[, c(1, 2)], "Time", "X1"); axis(1); axis(2)
make.plot(res[, c(3, 2)], "X2"); axis(1); axis(4)
title(main = "Van der Pol Oscillator", outer = TRUE)
par(oldpar)

