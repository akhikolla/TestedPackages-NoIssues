## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(comat)
data(raster_x, package = "comat")

## ---- echo=FALSE, fig.width=4, fig.height=4-----------------------------------
op = par(mar = rep(0, 4))
raster_x2 = apply(raster_x, 2, rev)
image(1:3, 1:3, t(raster_x2),
      col = c("#ffff64", "#006400", "#0046c8"),
      axes = FALSE, xlab = "", ylab = "")
par(op)

## ---- echo=FALSE--------------------------------------------------------------
raster_x

## -----------------------------------------------------------------------------
get_coma(raster_x)

## -----------------------------------------------------------------------------
get_coma(raster_x, neighbourhood = 8)

## -----------------------------------------------------------------------------
get_coma(raster_x, classes = c(1, 2))

## -----------------------------------------------------------------------------
my_coma = get_coma(raster_x)
get_cove(my_coma, normalization = "pdf")

