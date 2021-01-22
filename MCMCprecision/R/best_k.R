
#' Precision for the k Best-Performing Models
#'
#' Assesses the precision in estimating the ranking of the \eqn{k} best-performing models.
#'
#' @param samples a matrix with posterior samples (one per row) for the model posterior probabilities (one model per column). Can be estimated using \code{\link{stationary}} with the argument \code{summary = FALSE}.
#' @param k number of best-performing models to be considered
#' @param ties.method a character string specifying how ties are treated, see \code{\link[base]{rank}}
#'
#' @examples
#' # a sequence of uncorrelated model indices
#' mult <- rmultinom(1000, 1, c(.05, .6, .15, .12, .08))
#' idx <- apply(mult, 2, which.max)
#' z <- letters[idx]
#' stat <- stationary(z, summary = FALSE)
#' best_models(stat, 3)
#' @export
best_models <- function (samples,
                    k,
                    ties.method = "min"){

  if (k > ncol(samples)){
    warning ("'k' is larger then the number of models (i.e., the number of columns of 'samples').")
    k <- ncol(samples)
  }
  R <- nrow(samples)
  if (is.null(colnames(samples)))
    colnames(samples) <- paste0("M",1:ncol(samples))

  # 1. standard summary
  pp <- t(apply(samples, 2, summ.samples))

  # 2. select k best-performing models on average
  pp.sort <- sort(pp[,"50%"], decreasing = TRUE)
  best <- names(pp.sort)[1:k]
  overall <- rank(- pp[,"50%"], ties.method = ties.method)

  # 3. stability of rank order
  ranks <- apply(- samples, 1, rank, ties.method = ties.method)
  ranks.best <- ranks[best, , drop = FALSE]
  rank.stable <- apply(ranks.best == 1:k, 2, all)
  rank.identical <- rowMeans(ranks.best == 1:k)

  # 4. percentages: model rank <= k
  rank.mean <- rowMeans(ranks.best)
  rank.sd <- apply(ranks.best, 1, sd)
  in.best <- rowMeans(ranks.best <= k)

  # summarize
  summ <- cbind("#" = 1:k,
                "Prob." = pp[best,c(4,1,2)],
                "Rank.Mean" = rank.mean,
                "Rank.SD" = rank.sd,
                "Rank==#" = rank.identical,
                in.best)
  colnames(summ)[8] <- paste0("Rank<=",k)
  summ
}
