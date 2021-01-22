## ----setup, include = FALSE---------------------------------------------------
ev = Sys.getenv("USER") == "marius"
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = ev
)

## ----data_download, echo=TRUE, results='hide'---------------------------------
dest_dir = tempdir()
download.file("https://uni-muenster.sciebo.de/s/eP9E6OIkQbXrmsY/download", destfile=file.path(dest_dir, "MOD11A2.zip"),mode = "wb")
unzip(file.path(dest_dir, "MOD11A2.zip"), exdir = file.path(dest_dir,"MOD11A2"))
unlink(file.path(dest_dir, "MOD11A2.zip"))


## ---- eval=TRUE---------------------------------------------------------------
library(gdalcubes)
collection_formats()

## -----------------------------------------------------------------------------
img_col_file = tempfile(fileext=".db")
files = list.files(file.path(dest_dir,"MOD11A2"), pattern=".hdf$", full.names = TRUE)
create_image_collection(files, "MxD11A2", img_col_file)

## -----------------------------------------------------------------------------
x = raster_cube(image_collection(img_col_file))
x

bands(x)
dimensions(x)
srs(x)

## -----------------------------------------------------------------------------
srs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
MOD11A2.col = image_collection(img_col_file)
v = cube_view(srs=srs, extent=MOD11A2.col, nx = 400, dt="P1M", aggregation = "mean")
v
MOD11A2.cube = raster_cube(MOD11A2.col, v)
MOD11A2.cube

## -----------------------------------------------------------------------------
MOD11A2.bandselect = select_bands(MOD11A2.cube, c("LST_DAY","LST_NIGHT"))
MOD11A2.daynight_difference = apply_pixel(MOD11A2.bandselect, "0.02*(LST_DAY-LST_NIGHT)",names = "LST_difference")
MOD11A2.reduce = reduce_time(MOD11A2.daynight_difference, "median(LST_difference)")
plot(MOD11A2.reduce, col=heat.colors, key.pos=1)

## -----------------------------------------------------------------------------
library(magrittr) # get the pipe
v1 = cube_view(view=v, extent=list(t0="2018-06", t1="2018-09"))
raster_cube(image_collection(img_col_file), v1)  %>%
  select_bands(c("LST_DAY")) %>%
  apply_pixel("LST_DAY * 0.02") %>% # apply scale
  window_time(kernel=c(-1,1), window=c(1,0)) %>%
  plot(col=terrain.colors, key.pos=1, zlim=c(-25,25))

