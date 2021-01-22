### R code from vignette source 'stream_extension.Rnw'

###################################################
### code chunk number 1: stream_extension.Rnw:104-105
###################################################
options(width = 75, digits = 3, prompt = 'R> ', scipen = 3)


###################################################
### code chunk number 2: stream_extension.Rnw:197-198
###################################################
library("stream")


###################################################
### code chunk number 3: stream_extension.Rnw:201-216
###################################################
DSD_UniformNoise <- function(d = 2, range = NULL) {
  if(is.null(range)) range <- matrix(c(0, 1), ncol = 2, nrow = d,
    byrow = TRUE)
  structure(list(description = "Uniform Noise Data Stream", d = d,
    k = NA_integer_, range = range),
        class = c("DSD_UniformNoise", "DSD_R", "DSD"))
  }

get_points.DSD_UniformNoise <- function(x, n = 1,
  assignment = FALSE, ...) {
    data <- as.data.frame(t(replicate(n,
      runif(x$d, min = x$range[ , 1], max = x$range[ , 2]))))
    if(assignment) attr(data, "assignment") <- rep(NA_integer_, n)
    data
}


###################################################
### code chunk number 4: dsd_example
###################################################
stream <- DSD_UniformNoise()
stream
plot(stream, main = description(stream))


###################################################
### code chunk number 5: stream_extension.Rnw:246-248
###################################################
stream <- DSD_Gaussians(k = 1, d = 2, outliers = 1, space_limit = c(0,0.5),
                        outlier_options = list(outlier_horizon = 5))


###################################################
### code chunk number 6: stream_extension.Rnw:252-254
###################################################
points <- get_points(stream, n = 10, cluster = TRUE, outlier = TRUE)
points


###################################################
### code chunk number 7: stream_extension.Rnw:257-258
###################################################
attr(points, "cluster")


###################################################
### code chunk number 8: stream_extension.Rnw:261-262
###################################################
attr(points, "outlier")


###################################################
### code chunk number 9: stream_extension.Rnw:386-395 (eval = FALSE)
###################################################
## DSC_MyClusterer <- function(x) {
##   structure(
##     list(
##       description = "My new clusterer",
##       RObj = x
##     ), class = c("DSC_MyClusterer", "DSC_SinglePass", "DSC_Outlier",
##                  "DSC_Micro", "DSC_R", "DSC")
##   )
## }


###################################################
### code chunk number 10: stream_extension.Rnw:411-417 (eval = FALSE)
###################################################
## stream <- DSD_Gaussians(k = 1, d = 2, outliers = 1,
##                         space_limit = c(0, 1), variance_limit = .01,
##                         outlier_options = list(outlier_horizon = 20))
## points <- get_points(stream, n=20, cluster = TRUE, outlier = TRUE)
## dsc <- DSC_MyClusterer()
## assigns <- get_assignment(dsc, points, type="micro")


###################################################
### code chunk number 11: stream_extension.Rnw:419-420
###################################################
assigns <- readRDS("outlier_assignment.RDS")


