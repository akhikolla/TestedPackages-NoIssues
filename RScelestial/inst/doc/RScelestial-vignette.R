## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RScelestial)
# We load igraph for drawing trees. If you do not want to draw,
# there is no need to import igraph.
library(igraph)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("RScelestial")

## -----------------------------------------------------------------------------
# Following command generates ten samples with 20 loci. 
# Rate of mutations on each edge of the evolutionary tree is 1.5. 
D = synthesis(10, 20, 5, seed = 7)
D

## -----------------------------------------------------------------------------
seq = as.ten.state.matrix(D$seqeunce)
SP = scelestial(seq, return.graph = TRUE)
SP

## ----fig.width=5, fig.height=5------------------------------------------------
tree.plot(SP, vertex.size = 30)

## ----fig.width=5, fig.height=5------------------------------------------------
SP = scelestial(seq, root.assign.method = "fix", root = "C8", return.graph = TRUE)
tree.plot(SP, vertex.size = 30)

## ----fig.width=5, fig.height=5------------------------------------------------
SP = scelestial(seq, root.assign.method = "balance", return.graph = TRUE)
tree.plot(SP, vertex.size = 30)

## -----------------------------------------------------------------------------
D.distance.matrix <- distance.matrix.true.tree(D)
D.distance.matrix
SP.distance.matrix <- distance.matrix.scelestial(SP)
SP.distance.matrix
## Difference between normalized distance matrices
vertices <- rownames(SP.distance.matrix)
sum(abs(D.distance.matrix[vertices,vertices] - SP.distance.matrix))

