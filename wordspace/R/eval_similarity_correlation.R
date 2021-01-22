eval.similarity.correlation <- function (task, M, dist.fnc=pair.distances, details=FALSE, format=NA, taskname=NA, word1.name="word1", word2.name="word2", score.name="score", ...) {
  if (is.na(taskname)) taskname <- deparse(substitute(task))
  for (varname in c(word1.name, word2.name, score.name)) {
    if (!(varname %in% colnames(task))) stop(sprintf("gold standard does not have a column labelled '%s'", varname))
  }
  n.items <- nrow(task)

  w1 <- as.character(task[, word1.name])
  w2 <- as.character(task[, word2.name])
  if (!is.na(format)) {
    w1 <- convert.lemma(w1, format)
    w2 <- convert.lemma(w2, format)
  }
  score <- task[, score.name]
  distance <- dist.fnc(w1, w2, M, ...)   # vector of DSM distances between word pairs (can also be similarities)

  if (any(is.na(distance))) stop("missing values in distance/similarity vector must be replaced by +Inf or -Inf")
  is.pinf <- distance == Inf             # missing data points in case of distance 
  is.minf <- distance == -Inf            # missing data points in case of similarity
  is.missing <- is.pinf | is.minf
  n.missing <- sum(is.missing)
  if (all(is.missing)) stop("no distance/similarity values available, can't compute correlation")
  dist.range <- range(distance[!is.missing]) # fill in missing data with appropriate extreme of value range extended by 10%
  surrogates <- dist.range + diff(dist.range) * c(-.1, +.1) 
  distance[is.pinf] <- surrogates[2]
  distance[is.minf] <- surrogates[1]

  spearman <- cor.test(score, distance, method="spearman", exact=FALSE)
  rho <- abs(spearman$estimate)
  rho.pvalue <- spearman$p.value
  pearson <- cor.test(score, distance, method="pearson", conf.level=.95)
  r <- abs(pearson$estimate)
  r.confint <- if (pearson$estimate < 0) -rev(pearson$conf.int) else pearson$conf.int

  eval.result <- data.frame(
    rho = rho, p.value = rho.pvalue, missing = n.missing,
    r = r, r.lower = r.confint[1], r.upper = r.confint[2],
    row.names=taskname
  )

  if (details) {
    task$distance <- distance
    task$missing <- is.missing
    attr(task, "eval.result") <- eval.result
    attr(task, "taskname") <- taskname
    class(task) <- c("eval.similarity.correlation", "data.frame")
    task
  } else {
    eval.result
  }
}

print.eval.similarity.correlation <- function (x, ...) {
  NextMethod("print", x, ...) # print as data frame
  cat("Evaluation result:\n")
  print(attr(x, "eval.result"))
}

plot.eval.similarity.correlation <- function (x, y, line=TRUE, categories=NULL, cat.col=NA, cat.legend="bottomleft", pch=20, cex=1, xlim=NULL, ylim=NULL, xlab="human rating", ylab="distributional model", main=attr(x, "taskname"), ...) {
  if (!missing(categories)) {
    categories <- eval(substitute(categories), x, parent.frame())
    categories <- as.factor(categories)
    stopifnot(length(categories) == nrow(x))
    cat.types <- levels(categories)
    n.cat <- length(cat.types)
    if (is.na(cat.col)) {
      cat.col <- 1:n.cat
    } else {
      cat.col <- rep(cat.col, length.out=n.cat)
    }
    colours <- cat.col[categories]
  } else {
    n.cat <- 0
    colours <- rep("black", nrow(x))
  }
  
  if (is.null(xlim)) xlim <- range(x$score)
  if (is.null(ylim)) {
    ylim <- range(x$distance)
    ylim <- ylim + c(0, 0.05) * diff(ylim) # extend top of range by 5%
  }
  
  pch.vec <- ifelse(x$missing, 0, pch)
  cex.vec <- ifelse(x$missing, cex * 0.8, cex)
  plot(x$score, x$distance, col=colours, pch=pch.vec, cex=cex.vec, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)
  if (line) lines(lowess(x$score, x$distance), col="blue", lwd=3)

  if (n.cat > 0) {
    legend(cat.legend, inset=.02, bg="white", legend=cat.types, pch=pch, col=cat.col, pt.cex=cex*1.5)
  }
  
  viewport <- par("usr") # computer anchor coordinates in top left corner with 2% inset
  anchor.x <- sum(c(0.98, 0.02, 0, 0) * viewport)
  anchor.y <- sum(c(0, 0, 0.02, 0.98) * viewport)
  result <- attr(x, "eval.result")
  report <- with(result, sprintf("|rho| = %.3f, p = %.4f, |r| = %.3f .. %.3f", rho, p.value, r.lower, r.upper))
  if (result$missing > 0) report <- sprintf("%s (%d pairs not found)", report, result$missing)
  text(anchor.x, anchor.y, adj=c(0, 1), font=2, col="blue", labels=report)
}
