## determine majority class assignments for data points
majority.class <- function(clusters, gold) {
  class.dist <- table(clusters, gold) # class distribution for each cluster
  best.class <- colnames(class.dist)[ apply(class.dist, 1, which.max) ] # best class assignment for each cluster (by name)
  names(best.class) <- rownames(class.dist)
  best.class[ clusters ] # look up best class assignment for each data point
}

## clustering purity is easily computed from best assignments
cluster.purity <- function (clusters, gold) {
  stopifnot(length(gold) == length(clusters))
  best.class <- majority.class(clusters, gold)
  stopifnot(length(best.class) == length(clusters))
  100 * sum(as.character(best.class) == as.character(gold)) / length(clusters)
}

## clustering entropy (scaled to range 0..1 if scale=TRUE)
cluster.entropy <- function(clusters, gold, scale=FALSE) {
  stopifnot(length(gold) == length(clusters))
  class.dist <- table(clusters, gold)
  H.by.class <- apply(prop.table(class.dist, 1), 1,
                      function (p) sum(ifelse(p > 0, -p * log2(p), 0)))
  class.p <- apply(class.dist, 1, sum) / length(clusters)
  H <- sum(class.p * H.by.class) # weighted average of cluster entropies = mean clasification uncertainty
  if (scale) {
    p <- prop.table(table(gold)) # largest plausible clustering entropy = entropy of gold standard
    H <- H / sum(ifelse(p > 0, -p * log2(p), 0))
  }
  H
}

eval.clustering <- function (task, M, dist.fnc=pair.distances, ..., details=FALSE, format=NA, taskname=NA, scale.entropy=FALSE, n.clusters=NA, word.name="word", class.name="class") {
  if (is.na(taskname)) taskname <- deparse(substitute(task))
  if (!(word.name %in% colnames(task))) stop(sprintf("gold standard does not have a column labelled '%s'", word.name))
  if (!(class.name %in% colnames(task))) stop(sprintf("gold standard does not have a column labelled '%s'", class.name))
  
  orig.words <- as.character(task[, word.name])
  words <- if (is.na(format)) orig.words else convert.lemma(orig.words, format)
  gold <- as.factor(task[, class.name])
  n.gold <- length(unique(gold))
  if (is.na(n.clusters)) n.clusters <- n.gold

  ok <- words %in% rownames(M)          # check for missing words
  known.words <- words[ok]

  clusters <- rep("n/a", length(words))   # unknown words are assigned to a single cluster "n/a"
  names(clusters) <- words

  if (any(ok)) {
    ## create a dissimilarity structure for known words using pair.distances
    pairs <- t(combn(known.words, 2))
    distances <- dist.fnc(pairs[, 1], pairs[, 2], M, ...)
    class(distances) <- "dist"
    attr(distances, "Size") <- length(known.words)
    attr(distances, "Labels") <- known.words

    ## cluster assignments for known words, using the PAM algorithm
    known.clusters <- pam(distances, n.clusters, diss=TRUE, cluster.only=TRUE) 
    clusters[known.words] <- known.clusters # fill in cluster assignments for known words
  }
    
  best.class <- majority.class(clusters, gold) # best label for each cluster
  purity <- cluster.purity(clusters, gold)
  entropy <- cluster.entropy(clusters, gold, scale=scale.entropy)
  
  if (details) {
    res <- data.frame(
      word = orig.words,
      cluster = as.factor(clusters),
      label = factor(best.class, levels=levels(gold)),
      gold = gold,
      correct = best.class == gold,
      missing = !ok,
      stringsAsFactors = FALSE)
    attr(res, "purity") <- purity
    attr(res, "entropy") <- entropy
    attr(res, "taskname") <- taskname
    res
  } else {
    res <- data.frame(purity=purity, entropy=entropy, missing=sum(!ok))
    rownames(res) <- taskname
    res
  }
}
