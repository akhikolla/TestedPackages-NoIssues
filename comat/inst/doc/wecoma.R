## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(comat)
data(raster_x, package = "comat")
data(raster_w, package = "comat")

## -----------------------------------------------------------------------------
raster_x

## ---- echo=FALSE, fig.width=4, fig.height=4-----------------------------------
op = par(mar = rep(0, 4))
raster_x2 = apply(raster_x, 2, rev)
image(1:3, 1:3, t(raster_x2),
      col = c("#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", 
              "#E17C05", "#CC503E", "#94346E", "#6F4070", "#994E95", "#666666")[3:5],
      axes = FALSE, xlab = "", ylab = "")
par(op)

## -----------------------------------------------------------------------------
raster_w

## ---- echo=FALSE, fig.width=4, fig.height=4-----------------------------------
op = par(mar = rep(0, 4))
raster_w2 = apply(raster_w, 2, rev)
image(1:3, 1:3, t(raster_w2),
      col = rev(c("#3D1778", "#583A99", "#715DAA", "#8B7EBB", "#A69DCC",
              "#BFBADD", "#D8D4EC", "#EDECF9")),
      axes = FALSE, xlab = "", ylab = "")
par(op)

## -----------------------------------------------------------------------------
get_wecoma(raster_x, raster_w)

## -----------------------------------------------------------------------------
get_wecoma(raster_x, raster_w, fun = "focal", na_action = "omit")

## -----------------------------------------------------------------------------
my_wecoma = get_wecoma(raster_x, raster_w)
get_wecove(my_wecoma, normalization = "pdf")

