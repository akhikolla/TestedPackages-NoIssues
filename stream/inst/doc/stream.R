### R code from vignette source 'stream.Rnw'

###################################################
### code chunk number 1: stream.Rnw:138-139
###################################################
options(width = 75, digits = 3, prompt = 'R> ', scipen = 3)


###################################################
### code chunk number 2: stream.Rnw:682-685
###################################################
library("stream")
set.seed(1000)
stream <- DSD_Gaussians(k = 3, d = 2)


###################################################
### code chunk number 3: stream.Rnw:694-696
###################################################
dstream <- DSC_DStream(gridsize = .1, Cm = 1.2)
update(dstream, stream, n = 500)


###################################################
### code chunk number 4: initial_example
###################################################
km <- DSC_Kmeans(k = 3)
recluster(km, dstream)
plot(km, stream, type = "both")


###################################################
### code chunk number 5: stream.Rnw:912-917
###################################################
library("stream")
set.seed(1000)

stream <- DSD_Gaussians(k = 3, d = 3, noise = .05, p = c(.5, .3, .1))
stream


###################################################
### code chunk number 6: stream.Rnw:938-940
###################################################
p <- get_points(stream, n = 5)
p


###################################################
### code chunk number 7: stream.Rnw:951-953
###################################################
p <- get_points(stream, n = 100, class = TRUE)
head(p, n = 10)


###################################################
### code chunk number 8: static
###################################################
plot(stream, n = 500)


###################################################
### code chunk number 9: static_pc
###################################################
plot(stream, n = 500, method = "pc")


###################################################
### code chunk number 10: moa1
###################################################
set.seed(1000)
stream <- DSD_Benchmark(1)
stream


###################################################
### code chunk number 11: stream.Rnw:1011-1015 (eval = FALSE)
###################################################
## for(i in 1:4) {
##   plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
##   tmp <- get_points(stream, n = 1400)
## }


###################################################
### code chunk number 12: moa1
###################################################
plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
arrows(.15, .85, .85, .15, col = rgb(.8, .8, .8, .6), lwd = 10)
arrows(.15, .15, .85, .85, col = rgb(.8, .8, .8, .6), lwd = 10)
tmp <- get_points(stream, n = 1400)


###################################################
### code chunk number 13: moa2
###################################################
plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
arrows(.15, .85, .85, .15, col = rgb(.8, .8, .8, .6), lwd = 10)
arrows(.15, .15, .85, .85, col = rgb(.8, .8, .8, .6), lwd = 10)
tmp <- get_points(stream, n=1400)


###################################################
### code chunk number 14: moa3
###################################################
plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
arrows(.15,.85,.85,.15, col=rgb(.8,.8,.8,.6), lwd=10)
arrows(.15,.15,.85,.85, col=rgb(.8,.8,.8,.6), lwd=10)
tmp <- get_points(stream, n=1400)


###################################################
### code chunk number 15: moa4
###################################################
plot(stream, 250, xlim=c(0,1), ylim=c(0,1))
arrows(.15,.85,.85,.15, col=rgb(.8,.8,.8,.6), lwd=10)
arrows(.15,.15,.85,.85, col=rgb(.8,.8,.8,.6), lwd=10)


###################################################
### code chunk number 16: stream.Rnw:1068-1071 (eval = FALSE)
###################################################
## reset_stream(stream)
## animate_data(stream, n = 10000, horizon = 100,
##   xlim = c(0, 1), ylim = c(0, 1))


###################################################
### code chunk number 17: stream.Rnw:1077-1080 (eval = FALSE)
###################################################
## library("animation")
## animation::ani.options(interval = .1)
## ani.replay()


###################################################
### code chunk number 18: stream.Rnw:1087-1089 (eval = FALSE)
###################################################
## saveHTML(ani.replay())
## saveGIF(ani.replay())


###################################################
### code chunk number 19: stream.Rnw:1102-1107
###################################################
library("stream")
set.seed(1000)
stream <- DSD_Gaussians(k = 3, d = 2, outliers = 4,
                        outlier_options = list(outlier_horizon = 10000),
                        separation = 0.3, space_limit = c(0,1))


###################################################
### code chunk number 20: out1
###################################################
plot(stream, 10000, xlim = c(0, 1), ylim = c(0, 1))


###################################################
### code chunk number 21: stream.Rnw:1122-1125
###################################################
reset_stream(stream)
p <- get_points(stream, n = 10000, outlier = TRUE)
head(p)


###################################################
### code chunk number 22: stream.Rnw:1129-1131
###################################################
out_marks <- attr(p, "outlier")
sum(out_marks)


###################################################
### code chunk number 23: stream.Rnw:1134-1135
###################################################
which(out_marks)


###################################################
### code chunk number 24: stream.Rnw:1150-1156
###################################################
library("stream")
set.seed(1000)
stream1 <- DSD_Gaussians(k = 3, d = 2, variance_limit = 0.2,
                         space_limit = c(0, 5))
stream2 <- DSD_Gaussians(k = 3, d = 2, variance_limit = 2,
                         space_limit = c(0, 5))


###################################################
### code chunk number 25: dsd-lim1
###################################################
plot(stream1, 1000)


###################################################
### code chunk number 26: dsd-lim2
###################################################
plot(stream2, 1000)


###################################################
### code chunk number 27: stream.Rnw:1185-1191
###################################################
library("stream")
set.seed(1000)
stream1 <- DSD_Gaussians(k = 5, d = 2, variance_limit = 0.2,
                         space_limit = c(0, 7),
                         separation_type = "Mahalanobis",
                         separation = 4)


###################################################
### code chunk number 28: stream.Rnw:1193-1198
###################################################
set.seed(1000)
stream2 <- DSD_Gaussians(k = 5, d = 2, variance_limit = 0.2,
                         space_limit = c(0, 15),
                         separation_type = "Mahalanobis",
                         separation = 10)


###################################################
### code chunk number 29: dsd-ms1
###################################################
plot(stream1, 1000)


###################################################
### code chunk number 30: dsd-ms2
###################################################
plot(stream2, 1000)


###################################################
### code chunk number 31: stream.Rnw:1226-1233
###################################################
library("stream")
set.seed(1000)
stream1 <- DSD_Gaussians(k = 5, d = 2, outliers = 5, variance_limit = 0.2,
                         space_limit = c(0, 15), separation = 4,
                         separation_type = "Mahalanobis",
                         outlier_options = list(
                             outlier_virtual_variance = 0.3))


###################################################
### code chunk number 32: stream.Rnw:1235-1241
###################################################
set.seed(1000)
stream2 <- DSD_Gaussians(k = 5, d = 2, outliers = 5, variance_limit = 0.2,
                         space_limit = c(0, 40), separation = 4,
                         separation_type = "Mahalanobis",
                         outlier_options = list(
                             outlier_virtual_variance = 3))


###################################################
### code chunk number 33: dsd-mso1
###################################################
plot(stream1, 1000)


###################################################
### code chunk number 34: dsd-mso2
###################################################
plot(stream2, 1000)


###################################################
### code chunk number 35: stream.Rnw:1277-1280
###################################################
library("stream")
set.seed(1000)
stream <- DSD_Gaussians(k = 3, d = 5)


###################################################
### code chunk number 36: stream.Rnw:1285-1286 (eval = FALSE)
###################################################
## write_stream(stream, "data.csv", n = 100, sep = ",")


###################################################
### code chunk number 37: stream.Rnw:1320-1324
###################################################
file <- system.file("examples", "kddcup10000.data.gz", package = "stream")
stream_file <- DSD_ReadCSV(gzfile(file),
  take = c(1, 5, 6, 8:11, 13:20, 23:42), class = 42, k = 7)
stream_file


###################################################
### code chunk number 38: stream.Rnw:1337-1338
###################################################
get_points(stream_file, n = 5)


###################################################
### code chunk number 39: stream.Rnw:1345-1347
###################################################
stream_scaled <- DSD_ScaleStream(stream_file, center = TRUE, scale = TRUE)
get_points(stream_scaled, n = 5)


###################################################
### code chunk number 40: stream.Rnw:1378-1380
###################################################
data("EuStockMarkets", package = "datasets")
head(EuStockMarkets)


###################################################
### code chunk number 41: stream.Rnw:1387-1389
###################################################
replayer <- DSD_Memory(EuStockMarkets, k = NA)
replayer


###################################################
### code chunk number 42: stream.Rnw:1395-1397
###################################################
get_points(replayer, n = 5)
replayer


###################################################
### code chunk number 43: stream.Rnw:1405-1406 (eval = FALSE)
###################################################
## get_points(replayer, n = 2000)


###################################################
### code chunk number 44: stream.Rnw:1408-1410
###################################################
err <- try(get_points(replayer, n = 2000))
cat(err)


###################################################
### code chunk number 45: stream.Rnw:1424-1426
###################################################
reset_stream(replayer, pos = 100)
replayer


###################################################
### code chunk number 46: stream.Rnw:1804-1807
###################################################
library("stream")
set.seed(1000)
stream <- DSD_Gaussians(k = 3, d = 2, noise = .05)


###################################################
### code chunk number 47: stream.Rnw:1815-1817
###################################################
dstream <- DSC_DStream(gridsize = .1, Cm = 1.2)
dstream


###################################################
### code chunk number 48: stream.Rnw:1825-1827
###################################################
update(dstream, stream, n = 500)
dstream


###################################################
### code chunk number 49: stream.Rnw:1836-1837
###################################################
head(get_centers(dstream))


###################################################
### code chunk number 50: cluster
###################################################
plot(dstream, stream)


###################################################
### code chunk number 51: cluster-grid
###################################################
plot(dstream, stream, grid = TRUE)


###################################################
### code chunk number 52: stream.Rnw:2228-2232
###################################################
library("stream")
stream <- DSD_Gaussians(k = 3, d = 2, noise = .05)
dstream <- DSC_DStream(gridsize = .1)
update(dstream, stream, n = 2000)


###################################################
### code chunk number 53: stream.Rnw:2240-2241
###################################################
evaluate(dstream, stream, n = 100)


###################################################
### code chunk number 54: stream.Rnw:2252-2253
###################################################
evaluate(dstream, stream, measure = c("purity", "crand"), n = 500)


###################################################
### code chunk number 55: stream.Rnw:2294-2300
###################################################
set.seed(1000)
stream <- DSD_Benchmark(1)
dstream <- DSC_DStream(gridsize = .05, lambda = .01)
ev <- evaluate_cluster(dstream, stream,
  measure = c("numMicroClusters", "purity"), n = 5000, horizon = 100)
head(ev)


###################################################
### code chunk number 56: evaluation
###################################################
plot(ev[ , "points"], ev[ , "purity"], type = "l",
  ylab = "Avg. Purity", xlab = "Points")


###################################################
### code chunk number 57: stream.Rnw:2333-2338 (eval = FALSE)
###################################################
## set.seed(1000)
## stream <- DSD_Benchmark(1)
## dstream <- DSC_DStream(gridsize = .05, lambda = .01)
## r <- animate_cluster(dstream, stream, horizon = 100, n = 5000,
##      measure = "purity", plot.args = list(xlim = c(0, 1), ylim = c(0, 1)))


###################################################
### code chunk number 58: stream.Rnw:2368-2375
###################################################
library("stream")
set.seed(1000)
stream <- DSD_Gaussians(k = 3, d = 2, noise = .05)
dstream <- DSC_DStream(gridsize = .05, Cm = 1.5)

update(dstream, stream, n = 1000)
dstream


###################################################
### code chunk number 59: recluster
###################################################
plot(dstream, stream, type = "both")


###################################################
### code chunk number 60: recluster2
###################################################
km <- DSC_Kmeans(k = 3, weighted = TRUE)
recluster(km, dstream)
km
plot(km, stream, type = "both")


###################################################
### code chunk number 61: stream.Rnw:2436-2437
###################################################
evaluate(km, stream, measure = c("purity", "crand", "SSQ"), n = 1000)


###################################################
### code chunk number 62: stream.Rnw:2442-2444
###################################################
evaluate(km, stream, c(measure = "purity", "crand", "SSQ"), n = 1000,
  assign = "macro")


###################################################
### code chunk number 63: stream.Rnw:2467-2470
###################################################
points <- get_points(stream, n = 100)
assignment <- get_assignment(dstream, points, type = "macro")
assignment


###################################################
### code chunk number 64: silhouette
###################################################
assignment[is.na(assignment)] <- 0L
library("cluster")
plot(silhouette(assignment, dist = dist(points)))


###################################################
### code chunk number 65: stream.Rnw:2505-2517
###################################################
library("stream")
CustomCallback <- function() {
  env <- environment()
  all_measures <- c("LowestWeightPercentage")
  internal_measures <- c()
  external_measures <- all_measures
  outlier_measures <- c()
  this <- list(description = "Custom evaluation callback",
               env = environment())
  class(this) <- c("CustomCallback", "EvalCallback")
  this
}


###################################################
### code chunk number 66: stream.Rnw:2524-2534
###################################################
evaluate_callback.CustomCallback <- function(cb_obj, dsc, measure, points,
                                             actual, predict, outliers,
                                             predict_outliers,
                                             predict_outliers_corrid,
                                             centers, noise) {
    r <- list()
    if("LowestWeightPercentage" %in% measure)
        r$LowestWeightPercentage=min(get_weights(dsc))/sum(get_weights(dsc))
    r
}


###################################################
### code chunk number 67: stream.Rnw:2545-2551
###################################################
stream <- DSD_Gaussians(k = 3, d = 2, p = c(0.2, 0.4, 0.4))
km <- DSC_Kmeans(3)
update(km, stream, n=500)
evaluate_with_callbacks(km, stream, type="macro", n=500,
                        measure = c("crand","LowestWeightPercentage"),
                        callbacks = list(cc=CustomCallback()))


###################################################
### code chunk number 68: data_bng
###################################################
set.seed(1000)
library("stream")
stream <- DSD_Memory(DSD_BarsAndGaussians(noise = .05), n = 1500)
stream
plot(stream)


###################################################
### code chunk number 69: stream.Rnw:2745-2753
###################################################
algorithms <- list(
  'Sample' = DSC_TwoStage(micro = DSC_Sample(k = 100),
    macro = DSC_Kmeans(k = 4)),
  'Window' = DSC_TwoStage(micro = DSC_Window(horizon = 100),
    macro = DSC_Kmeans(k = 4)),
  'D-Stream' = DSC_DStream(gridsize = .7, Cm = 1.5),
  'DBSTREAM' = DSC_DBSTREAM(r = .45)
)


###################################################
### code chunk number 70: stream.Rnw:2763-2767
###################################################
for(a in algorithms) {
  reset_stream(stream)
  update(a, stream, n = 1000)
}


###################################################
### code chunk number 71: stream.Rnw:2772-2773
###################################################
sapply(algorithms, nclusters, type = "micro")


###################################################
### code chunk number 72: microclusters
###################################################
op <- par(no.readonly = TRUE)
layout(mat = matrix(1:length(algorithms), ncol = 2))
for(a in algorithms) {
  reset_stream(stream)
  plot(a, stream, main = description(a), type = "micro")
}
par(op)


###################################################
### code chunk number 73: microclusters_assignment
###################################################
op <- par(no.readonly = TRUE)
layout(mat = matrix(1:length(algorithms), ncol = 2))
for(a in algorithms) {
  reset_stream(stream)
  plot(a, stream, main = description(a),
    assignment = TRUE, weight = FALSE, type = "micro")
}
par(op)


###################################################
### code chunk number 74: stream.Rnw:2851-2858
###################################################
sapply(algorithms, FUN=function(a) {
  reset_stream(stream, pos = 1001)
  evaluate(a, stream,
    measure = c("numMicroClusters", "purity"),
    type = "micro",
    n = 500)
})


###################################################
### code chunk number 75: macroclusters
###################################################
op <- par(no.readonly = TRUE)
layout(mat=matrix(1:length(algorithms), ncol = 2))
for(a in algorithms) {
  reset_stream(stream)
  plot(a, stream, main = description(a), type = "both")
}
par(op)


###################################################
### code chunk number 76: stream.Rnw:2902-2908
###################################################
sapply(algorithms, FUN = function(a) {
  reset_stream(stream, pos = 1001)
  evaluate(a, stream, measure = c("numMacroClusters", "purity",
      "SSQ", "cRand", "silhouette"),
    n = 500, assign = "micro", type = "macro")
})


###################################################
### code chunk number 77: stream.Rnw:2930-2932
###################################################
set.seed(0)
stream <- DSD_Memory(DSD_Benchmark(1), n = 5000)


###################################################
### code chunk number 78: stream.Rnw:2942-2950
###################################################
algorithms <- list(
  'Sample' = DSC_TwoStage(micro = DSC_Sample(k = 100, biased = TRUE),
    macro = DSC_Kmeans(k = 2)),
  'Window' = DSC_TwoStage(micro = DSC_Window(horizon = 100, lambda = .01),
    macro = DSC_Kmeans(k = 2)),
  'D-Stream' = DSC_DStream(gridsize = .1, lambda = .01),
  'DBSTREAM' = DSC_DBSTREAM(r = .05, lambda = .01)
)


###################################################
### code chunk number 79: stream.Rnw:2961-2966
###################################################
evaluation <- lapply(algorithms, FUN = function(a) {
  reset_stream(stream)
  evaluate_cluster(a, stream, horizon = 100, n = 5000, measure = "crand",
    type = "macro", assign = "micro")
})


###################################################
### code chunk number 80: stream.Rnw:2982-2985
###################################################
Position <- evaluation[[1]][ , "points"]
cRand <- sapply(evaluation, FUN = function(x) x[ , "cRand"])
head(cRand)


###################################################
### code chunk number 81: dynamic
###################################################
matplot(Position, cRand, type = "l", lwd = 2)
legend("bottomleft", legend = names(evaluation),
  col = 1:6, lty = 1:6, bty = "n", lwd = 2)


###################################################
### code chunk number 82: dynamic_box
###################################################
boxplot(cRand, las = 2, cex.axis = .8)


###################################################
### code chunk number 83: stream.Rnw:3043-3051 (eval = FALSE)
###################################################
## library("stream")
## con <- gzcon(
##   url(paste("http://archive.ics.uci.edu/ml/machine-learning-databases/",
##     "kddcup99-mld/kddcup.data.gz", sep="")))
## 
## stream <- DSD_ReadCSV(con, take=c(1, 5, 6, 8:11, 13:20, 23:42),
##     class=42, k=7)
## stream2 <- DSD_ScaleStream(stream, n=1000)


###################################################
### code chunk number 84: stream.Rnw:3057-3059 (eval = FALSE)
###################################################
## dstream <- DSC_DStream(gridsize = .5, gaptime = 10000L, lambda = .01)
## update(dstream, stream2, n = 4000000, verbose = TRUE)


