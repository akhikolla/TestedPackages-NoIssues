plot.dist.matrix <- function (x, y, labels=rownames(x), show.labels=TRUE, label.pos=3, selected=attr(x, "selected"), show.selected=TRUE, col="black", cex=1, pch=20, pt.cex=1.2, selected.cex=1.2, selected.col="red", show.edges=TRUE, edges.lwd=6, edges.col="#AABBFF", edges.threshold=quantile(x, 2/3), method=c("isomds", "sammon"), aspect=1, expand=.05, ...) {
  stopifnot(inherits(x, "dist.matrix"))
  if (!missing(y)) stop("plot.dist.matrix() doesn't take a second argument (y)")
  if (!isTRUE(attr(x, "symmetric"))) stop("only symmetric distance matrices can be plotted")
  if (isTRUE(attr(x, "similarity"))) stop("similarity matrices are not supported. Please provide a distance matrix for this plot.")
  method <- match.arg(method)
  if (is.null(labels)) {
    show.labels <- FALSE
  } else {
    if (length(labels) != nrow(x)) stop("wrong number of labels specified")
  }
    
  coords <- if (method == "isomds") isoMDS(x, k=2, trace=FALSE)$points else sammon(x, k=2, trace=FALSE)$points
  x.range <- extendrange(coords[, 1], f=expand)
  y.range <- extendrange(coords[, 2], f=expand)
  .asp <- diff(x.range) / diff(y.range)
  if (.asp < aspect) {
    x.range <- extendrange(x.range, f=(aspect/.asp-1)/2)
  } else if (.asp > aspect) {
    y.range <- extendrange(y.range, f=(.asp/aspect-1)/2)
  }

  ## set up plot region
  plot(coords, type="n", xlim=x.range, ylim=y.range, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ...)

  ## draw edges if requested
  if (show.edges) {
    midx <- t(combn(1:nrow(x), 2))   # 2d-index for upper triangle of distance matrix
    P1 <- midx[, 1]                  # idx of start node of edge
    P2 <- midx[, 2]                  # idx of end node of edge
    len <- x[midx]                   # distance between P1 and P2
    lwd.vec <- edges.lwd * (1 - len / edges.threshold)
    idx <- lwd.vec > 0.1             # draw only edges with minimum lwd of 0.1
    segments(coords[P1[idx], 1], coords[P1[idx], 2], coords[P2[idx], 1], coords[P2[idx], 2], lwd=lwd.vec[idx], col=edges.col)
  }

  ## draw points
  if (is.null(selected) || !show.selected) selected <- rep(FALSE, nrow(x))
  col.vec <- ifelse(selected, selected.col, col)
  cex.vec <- cex * pt.cex * ifelse(selected, selected.cex, 1)
  points(coords, pch=pch, col=col.vec, cex=cex.vec)

  ## draw labels if requested
  if (show.labels) {
    cex.vec <- cex * ifelse(selected, selected.cex, 1)
    text(coords[, 1], coords[, 2], labels=labels, pos=label.pos, font=2, cex=cex.vec, col=col.vec)
  }

  ## return matrix of MDS coordinates with row labels
  if (!is.null(labels)) rownames(coords) <- labels
  invisible(coords)
}
