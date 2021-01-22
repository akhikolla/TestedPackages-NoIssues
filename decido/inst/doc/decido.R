## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----example------------------------------------------------------------------
library(decido)
x <- c(0, 0, 0.75, 1, 0.5, 0.8, 0.69)
y <- c(0, 1, 1, 0.8, 0.7, 0.6, 0)
(ind <- earcut(cbind(x, y)))


plot_ears(cbind(x, y), ind)

## -----------------------------------------------------------------------------
## polygon with a hole
x <- c(0, 0, 0.75, 1, 0.5, 0.8, 0.69,
     0.2, 0.5, 0.5, 0.3, 0.2)
y <- c(0, 1, 1, 0.8, 0.7, 0.6, 0,
     0.2, 0.2, 0.4, 0.6, 0.4)
ind <- earcut(cbind(x, y), holes = 8)
plot_ears(cbind(x, y), ind)

## -----------------------------------------------------------------------------
plot_ears(cbind(x, y), ind, col = "grey", border = NA)
text(x, y, labels = seq_along(x), pos = 2)


## -----------------------------------------------------------------------------
## add another hole
x <- c(0, 0, 0.75, 1, 0.5, 0.8, 0.69,
     0.2, 0.5, 0.5, 0.3, 0.2,
      0.15, 0.23, 0.2)
y <- c(0, 1, 1, 0.8, 0.7, 0.6, 0,
     0.2, 0.2, 0.4, 0.6, 0.4,
      0.65, 0.65, 0.81)
ind <- earcut(cbind(x, y), holes = c(8, 13))
plot_ears(cbind(x, y), ind, col = "grey")


## ----headers, eval = FALSE, include = TRUE------------------------------------
#  ## this code is not run in the vignette
#  ## but can be run as is with decido installed
#  library(Rcpp)
#  
#  cppFunction(
#    depends = "decido"
#    , includes = '#include "decido/decido.hpp"'
#    , code = '
#      Rcpp::IntegerVector earcut0( SEXP polygon ) {
#        return decido::api::earcut( polygon );
#      }
#    '
#  )
#  
#  poly <- list(matrix(c(0,0,0,1,1,1,1,0,0,0), ncol = 2, byrow = T))
#  earcut0( poly )
#  # [1] 1 4 3 3 2 1

## ----triplot------------------------------------------------------------------
## plot triangles (input is triplets of xy coordinates)
## with one side an oriented arrow
plot_tri <- function(x, add = TRUE) {
  if (!add) plot(x, asp = 1, type = "n")
  idx <- c(1:3, 1)
  for (i in seq_len(nrow(x)/3)) {
    #print(idx)
    ## plot only the first arrow
    arrows(x[idx[1], 1], x[idx[1], 2], 
           x[idx[2], 1], x[idx[2], 2], length = 0.25, lwd = 2, col = "firebrick")
    ## segments the rest
    segments(x[idx[2:3], 1], x[idx[2:3], 2], 
           x[idx[3:4], 1], x[idx[3:4], 2])
  
    idx <- idx + 3
    
  }
}


## ----orientation--------------------------------------------------------------
plot_tri(cbind(x, y)[ind, ], add = FALSE)

## -----------------------------------------------------------------------------
x <- c(0, 0, 1, 1,
       0.4, 0.2, 0.2, 0.4,
       0.6, 0.8, 0.8, 0.6
)
y <- c(0, 1, 1, 0,
       0.2, 0.2, 0.4, 0.4,
       0.6, 0.6, 0.4, 0.4
)
ind <- decido::earcut(cbind(x, y), holes = c(5, 9))
plot_ears(cbind(x, y), ind, col = "grey")
title("triangle plot, two holes")
plot_holes(cbind(x, y), holes = c(5, 9), col = "grey")
title("path plot, two holes")

ind <- decido::earcut(cbind(x, y), holes = 5)
plot_ears(cbind(x, y), ind, col = "grey")
title("triangle plot, two holes as one")
plot_holes(cbind(x, y), holes = 5, col = "grey")
title("path plot, two holes as one")

## -----------------------------------------------------------------------------
library(oz)
oz_ring <- oz::ozRegion(states = FALSE)
ring <- oz_ring$lines[[6]]
indices <- earcut(ring[c("x", "y")])
plot_ears(cbind(ring$x, ring$y), indices)

## -----------------------------------------------------------------------------
vecrot <- function(x, k) {
  if (k < 0 || k > length(x)) warning("k out of bounds of 'x' index")
  k <- k %% length(x)
  ## from DescTools::VecRot
  rep(x, times = 2)[(length(x) - k + 1):(2 * length(x) - k)]
}

## -----------------------------------------------------------------------------
x <- c(0, 0, 0.75, 1, 0.5, 0.8, 0.69)
y <- c(0, 1, 1, 0.8, 0.7, 0.6, 0)
(ind <- earcut(cbind(x, y)))

plot_ears(cbind(x, y), ind)
points(x[1], y[1], pch = 19, col = "firebrick")
title("original")
op <- par(mfrow = c(2, 3))
for (rot in head(seq_along(x), -1)) {
xx <- vecrot(x, rot); yy <- vecrot(y, rot)
ind <- earcut(cbind(xx, yy))
plot_ears(cbind(xx, yy), ind)
title(sprintf("rot %i", rot))
points(xx[1], yy[1], pch = 19, col = "firebrick")
}
par(op)



## ---- include=FALSE-----------------------------------------------------------
# , and it's also not the only way to do it, the grid package uses a grouping vector rather than a sparse index like this. The spatial packages sp and sf explicitly use structural hierarchies rather than more abstract specifications. The fortify-approach in ggplot2 is more like the grid one. A sparse representation is closer to what is needed for topological operations and visualization, consider that when we have triangles there are no need for "holes", we can identify which triangles will be plotted and how (or not), and an index into the vertices available becomes a key efficiency feature. (See [silicate](https://github.com/hypertidy/silicate) for a lot more on this topic). 

