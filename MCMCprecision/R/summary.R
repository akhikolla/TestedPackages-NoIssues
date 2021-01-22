#' Summary for Posterior Model Probabilities
#'
#' Summary for a sample of posterior model probabilities (\code{\link{stationary}}).
#' Also provides the estimated effective sample size and summaries for all pairwise Bayes factors.
#'
#' @param object posterior samples of the stationary distribution (rows = samples; columns = models).
#' @param BF whether to compute summaries for all pairwise Bayes factors.
#' @param logBF whether to summarize log Bayes factors instead of Bayes factors.
#' @param ... passed to \code{\link{fit_dirichlet}} to estimate effective sample size
#'     (e.g., \code{maxit} and \code{abstol}).
#'
#' @details
#' Effective sample is estimated by fitting a Dirichlet model to the
#' posterior model probabilities (thereby assuming that samples were drawn from
#' an equivalent multinomial distribution based on independent sampling).
#' More precisely, sample size is estimated by the sum of the Dirichlet parameters
#' \eqn{\sum\alpha[i]} minus the prior sample size \eqn{\epsilon*M^2}
#' (where \eqn{M} is the number of sampled models and \eqn{\epsilon} the
#' prior parameter for each transition frequency).
#'
#' @examples
#' P <- matrix(c(.1,.5,.4,
#'               0, .5,.5,
#'               .9,.1,0), ncol = 3, byrow=TRUE)
#' z <- rmarkov(1000, P)
#' samples <- stationary(z, summary = FALSE)
#' summary(samples)
#'
#' @return a list with estimates for
#'     \code{"pp"} = posterior model probabilities,
#'     \code{"n.eff"} = effective sample size, and
#'     \code{"bf"} = pairwise Bayes factors (optional)
#'
#' @seealso \code{\link{stationary}}, \code{\link{fit_dirichlet}}
#' @export
summary.stationary <- function(object, BF = FALSE, logBF = FALSE, ...){
  samples <- object
  pp <- t(apply(samples, 2, summ.samples))

  n.eff <- NA
  if (attr(samples, "method") != "iid"){
    try ({
      est <- fit_dirichlet(na.omit(samples), ...)$sum
      n.eff <-  floor(est - attr(object, "epsilon")*ncol(samples)^2)
    })
  }

  res <- list(pp=pp, n.eff=n.eff)

  if (BF){
    combs <- combn(ncol(samples),2)

    if (!is.null(colnames(samples))){
      rownames(pp) <- colnames(samples)
      c1 <- rownames(pp)[combs[1,]]
      c2 <- rownames(pp)[combs[2,]]
      bf.names <- paste0("BF_", apply(rbind(c1,c2),2,paste,collapse=""))

    } else {
      rownames(pp) <- paste0("M",1:ncol(samples))
      bf.names <- paste0("BF_", apply(combs,2,paste,collapse=""))
    }

    res$bf <- matrix(NA, ncol(combs), 5,
                     dimnames = list(BF = bf.names,
                                     Statistic = colnames(pp)))
    for (bb in 1:ncol(combs)){
      suppressWarnings(
        bf.tmp <- log(samples[,combs[1,bb]])-log(samples[,combs[2,bb]])
      )
      if (!logBF)
        bf.tmp <- exp(bf.tmp)
      res$bf[bb,] <- summ.samples(bf.tmp)
    }
  }
  return (res)
}

summ.samples <- function(x, probs = c(.05,.5,.95)){
  c(Mean=mean(x, na.rm = TRUE),
    SD=sd(x, na.rm = TRUE),
    quantile(x, probs = probs, na.rm = TRUE))
}
