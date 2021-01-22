eval.multiple.choice <-
function (task, M, dist.fnc=pair.distances, ..., details=FALSE, format=NA, taskname=NA, target.name="target", correct.name="correct", distractor.name="^distract") {
  if (is.na(taskname)) taskname <- deparse(substitute(task))
  if (!(target.name %in% colnames(task))) stop(sprintf("gold standard does not have a column labelled '%s'", target.name))
  if (!(correct.name %in% colnames(task))) stop(sprintf("gold standard does not have a column labelled '%s'", correct.name))
  idx.distract <- grepl(distractor.name, colnames(task), perl=TRUE)
  if (!any(idx.distract)) stop(sprintf("no distractors matching /%s/ found in gold standard", distractor.name))

  targets <- as.character(task[, target.name])
  choices <- cbind(as.character(task[, correct.name]),
                   as.matrix(task[, idx.distract, drop=FALSE]))
  mode(choices) <- "character" # first column contains correct choice, further columns are distractors
  n.choices <- ncol(choices)
  n.items <- nrow(task)
  
  w1 <- rep(targets, n.choices)
  w2 <- as.vector(choices)
  if (!is.na(format)) {
    w1 <- convert.lemma(w1, format)
    w2 <- convert.lemma(w2, format)
  }
  distance <- dist.fnc(w1, w2, M, ...)
  is.similarity <- isTRUE(attr(distance, "similarity"))
  distance <- matrix(distance, ncol=n.choices, byrow=FALSE) # distances in matrix format corresponding to choices
  if (is.similarity) distance <- -distance # simple trick 
  
  res.list <- lapply(1:n.items, function (i) {
    d <- distance[i, ]
    ranks <- rank(d, ties.method="max") # so we don't get a correct answer if it is tied with a distractor
    best <- which.min(ranks)
    correct <- if (d[best] < Inf) ranks[1] == 1 else NA # whether correct answer is ranked first (NA if all pairings not in DSM)
    data.frame(
      target=targets[i], correct=correct,
      best.choice=choices[i, best], best.dist=d[best],
      correct.choice=choices[i, 1], correct.rank=ranks[1], correct.dist=d[1],
      row.names=NULL, stringsAsFactors=FALSE
    )
  })
  res <- do.call("rbind", res.list)
  
  if (details) {
    if (is.similarity) {
      res$best.dist    <- -res$best.dist    # convert back to similarity scores
      res$correct.dist <- -res$correct.dist # (can't use transform() because of warnings from "make check")
    }
    res
  } else {
    tp <- sum(res$correct, na.rm=TRUE)
    data.frame(
      accuracy=100 * tp / n.items,
      TP=tp, FP=n.items - tp, missing=sum(distance[, 1] == Inf),
      row.names=taskname
    )
  }
}
