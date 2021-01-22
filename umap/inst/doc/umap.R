## ---- echo=FALSE--------------------------------------------------------------
## block with some startup/background objects functions
library(umap)

plot.iris = function(x, labels,
         main="A UMAP visualization of the Iris dataset",
         colors=c("#ff7f00", "#e377c2", "#17becf"),
         pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=0.85) {

  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)

  labels.u = unique(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }

  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

set.seed(123456)

## -----------------------------------------------------------------------------
head(iris, 3)

## -----------------------------------------------------------------------------
iris.data = iris[, grep("Sepal|Petal", colnames(iris))]
iris.labels = iris[, "Species"]

## ----iris.umap----------------------------------------------------------------
library(umap)
iris.umap = umap(iris.data)

## ----umap.print---------------------------------------------------------------
iris.umap

## ----umap.layout--------------------------------------------------------------
head(iris.umap$layout, 3)

## ---- fig.width=3.2, fig.height=3.2, dpi=150----------------------------------
plot.iris(iris.umap, iris.labels)

## -----------------------------------------------------------------------------
iris.wnoise = iris.data + matrix(rnorm(150*40, 0, 0.1), ncol=4)
colnames(iris.wnoise) = colnames(iris.data)
head(iris.wnoise, 3)

## -----------------------------------------------------------------------------
iris.wnoise.umap = predict(iris.umap, iris.wnoise)
head(iris.wnoise.umap, 3)

## ---- fig.width=3.6, fig.height=3.6, dpi=150----------------------------------
plot.iris(iris.umap, iris.labels)
plot.iris(iris.wnoise.umap, iris.labels, add=T, pch=4, legend.suffix=" (with noise)")

## ----defaults, eval=FALSE-----------------------------------------------------
#  umap.defaults

## ----defaults2, eval=TRUE, echo=FALSE, collapse=TRUE--------------------------
umap.defaults

## ----custom.config, eval=TRUE-------------------------------------------------
custom.config = umap.defaults
custom.config$random_state = 123

## ----custom2, fig.width=3.6, fig.height=3.6, dpi=150--------------------------
iris.umap.config = umap(iris.data, config=custom.config)
plot.iris(iris.umap.config, iris.labels,
          main="Another UMAP visualization (different seed)")

## ----custom3, eval=FALSE------------------------------------------------------
#  iris.umap.args = umap(iris.data, random_state=123)

## -----------------------------------------------------------------------------
iris.dist = as.matrix(dist(iris.data))
iris.umap.dist = umap(iris.dist, config=custom.config, input="dist")
iris.umap.dist

## -----------------------------------------------------------------------------
iris.umap$knn

## -----------------------------------------------------------------------------
# extract information on 10 nearest neighbors from iris.umap
iris.neighbors = iris.umap$knn$indexes[, 1:10]
iris.neighbors.distances = iris.umap$knn$distances[, 1:10]
# construct an object with indexes and distances
iris.knn.10 = umap.knn(indexes=iris.neighbors, distances=iris.neighbors.distances)
iris.knn.10
# perform an embedding using the precomputed nearest neighbors
iris.umap.knn = umap(iris.data, config=custom.config, n_neighbors=10, knn=iris.knn.10)

## ---- eval=TRUE---------------------------------------------------------------
# predict in batch, display first item
predict(iris.umap, iris.wnoise)[1, , drop=FALSE]
# predict only first item
predict(iris.umap, iris.wnoise[1,,drop=FALSE])

## ---- eval=FALSE--------------------------------------------------------------
#  iris.umap.learn = umap(iris.data, method="umap-learn")

## ----show.plot.iris-----------------------------------------------------------
plot.iris

## -----------------------------------------------------------------------------
sessionInfo()

