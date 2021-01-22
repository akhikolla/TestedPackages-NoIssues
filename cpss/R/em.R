#' Detecting changes in exponential family
#'
#' @param dataset a numeric matrix of dimension \eqn{n\times d}, where each row represents an observation and each column stands for a variable. A numeric vector could also be acceptable for univariate observations.
#' @param family a character string indicating the underlying distribution. Currently, detecting changes in binomial ("\code{binom}"), multinomial ("\code{multinom}"), Poisson ("\code{pois}"), exponential ("\code{exp}"), geometric ("\code{geom}"), dirichlet ("\code{diri}"), gamma ("\code{gamma}"), beta ("\code{beta}"), chi-square ("\code{chisq}") and inverse gaussian ("\code{invgauss}") distributions are supported.
#' @param size an integer indicating the number of trials if \code{family} = "binom" or \code{family} = "multinom".
#' @param algorithm a character string specifying the change-point searching algorithm, one of four state-of-the-art candidates "SN" (segment neighborhood), "BS" (binary segmentation), "WBS" (wild binary segmentation) and "PELT" (pruned exact linear time) algorithms.
#' @param dist_min an integer indicating the minimum distance between two successive candidate change-points, with a default value \eqn{floor(log(n))}.
#' @param ncps_max an integer indicating the maximum number of change-points searched for, with a default value \eqn{ceiling(n^0.4)}.
#' @param pelt_pen_val a numeric vector specifying the collection of candidate values of the penalty if the "PELT" algorithm is used.
#' @param pelt_K a numeric value to adjust the pruning tactic, usually is taken to be 0 if negative log-likelihood is used as a cost; more details can be found in Killick et al. (2012).
#' @param wbs_nintervals an integer indicating the number of random intervals drawn in the "WBS" algorithm and a default value 500 is used.
#' @param criterion a character string indicating which model selection criterion, "cross- validation" ("CV") or "multiple-splitting" ("MS"), is used.
#' @param times an integer indicating how many times of sample-splitting should be performed; if "CV" criterion is used, it should be set as 2.
#'
#' @return \code{cpss.em} returns an object of an \proglang{S4} class, called "\code{cpss}", which collects data and information required for further change-point analyses and summaries. See \code{\link{cpss.custom}}.
#' @export
#'
#' @references
#' Killick, R., Fearnhead, P., and Eckley, I. A. (2012). Optimal Detection of Changepoints With a Linear Computational Cost. Journal of the American Statistical Association, 107(500):1590–1598.
#' @seealso \code{\link{cpss.meanvar}} \code{\link{cpss.mean}} \code{\link{cpss.var}}
#' @examples
#' library("cpss")
#' set.seed(666)
#' n <- 1000
#' tau <- c(100, 300, 700, 900)
#' tau_ext <- c(0, tau, n)
#' theta <- c(1, 0.2, 1, 0.2, 1)
#' seg_len <- diff(c(0, tau, n))
#' y <- unlist(lapply(seq(1, length(tau) + 1), function(k) {
#'   rexp(seg_len[k], theta[k])
#' }))
#' res <- cpss.em(
#'   y, family = "exp", algorithm = "WBS",
#'   dist_min = 10, ncps_max = 10,
#'   criterion = "MS", times = 10
#' )
#' cps(res)
#' # [1] 100 299 705 901
cpss.em <- function(dataset, family, size = NULL, algorithm = "BS", dist_min = floor(log(n)), ncps_max = ceiling(n^0.4), pelt_pen_val = NULL, pelt_K = 0, wbs_nintervals = 500, criterion = "CV", times = 2) {

  dataset <- as.matrix(dataset)
  n <- nrow(dataset)
  param.opt <- list()
  param.opt$em <- family
  param.opt$N <- size
  if (family %in% c("binom", "multinom", "pois", "exp", "geom")) {
    res <- cpss.custom(
      dataset,
      n,
      g_subdat,
      g_param_em,
      g_cost_em,
      algorithm,
      dist_min,
      ncps_max,
      pelt_pen_val,
      pelt_K,
      wbs_nintervals,
      criterion,
      times,
      model = "em",
      g_smry = g_smry_em,
      easy_cost = NULL,
      param.opt = param.opt
    )
  } else if (family %in% c("diri", "gamma", "beta", "chisq", "invgauss")) {
    res <- cpss.custom(
      dataset,
      n,
      g_subdat,
      g_param_em,
      g_cost_em,
      algorithm,
      dist_min,
      ncps_max,
      pelt_pen_val,
      pelt_K,
      wbs_nintervals,
      criterion,
      times,
      model = "em",
      g_smry = NULL,
      easy_cost = NULL,
      param.opt = param.opt
    )
  } else {
    stop("Not supported member of exponential family!")
  }

  res@call <- list(call = match.call())

  return(res)
}
