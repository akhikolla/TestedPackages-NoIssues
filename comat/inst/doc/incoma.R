## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(comat)
data(raster_x, package = "comat")
data(raster_y, package = "comat")

## -----------------------------------------------------------------------------
raster_x

## ---- echo=FALSE, fig.width=4, fig.height=4-----------------------------------
op = par(mar = rep(0, 4))
raster_x2 = apply(raster_x, 2, rev)
image(1:3, 1:3, t(raster_x2),
      col = c("#ffff64", "#006400", "#BE9600"),
      axes = FALSE, xlab = "", ylab = "")
par(op)

## -----------------------------------------------------------------------------
raster_y

## ---- echo=FALSE, fig.width=4, fig.height=4-----------------------------------
op = par(mar = rep(0, 4))
raster_y2 = apply(raster_y, 2, rev)
image(1:3, 1:3, t(raster_y2),
      col = c("#fcffc2", "#724501"),
      axes = FALSE, xlab = "", ylab = "")
par(op)

## -----------------------------------------------------------------------------
get_incoma(list(raster_x, raster_y))

## -----------------------------------------------------------------------------
my_incoma = get_incoma(list(raster_x, raster_y))
get_wecove(my_incoma, normalization = "pdf")

