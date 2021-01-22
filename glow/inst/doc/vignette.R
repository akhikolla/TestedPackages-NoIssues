## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("traversc/glow")

## ----eval=FALSE---------------------------------------------------------------
#  library(glow)
#  library(ggplot2)
#  library(viridisLite) # Magma color scale
#  
#  # Load the dataset
#  data(diamonds)
#  gm <- GlowMapper$new(xdim=2000, ydim = 2000, blend_mode = "screen", nthreads=nt)
#  gm$map(x=diamonds$carat, y=diamonds$price, intensity=1, radius = .1)
#  pd <- gm$output_dataframe(saturation = 1)
#  ggplot() +
#    geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
#    scale_fill_gradientn(colors = additive_alpha(magma(12))) +
#    coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) +
#    labs(x = "carat", y = "price") +
#    theme_night(bgcolor = magma(1))

