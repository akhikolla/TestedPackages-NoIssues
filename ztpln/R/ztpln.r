#' The zero-truncated compund poisson-lognormal distributions and their
#' mixtures
#'
#' Functions for obtaining the density, random deviates and maximum likelihood
#' estimates of the zero-truncated Poisson lognormal distributions and their
#' mixtures.
#'
#' @author Masatoshi Katabuchi <mattocci27@gmail.com> 
#'
#' @references Bulmer, M. G. 1974. On Fitting the Poisson Lognormal Distribution to Species-Abundance Data. Biometrics 30:101-110.
#' @references Inouye, D., E. Yang, G. Allen, and P. Ravikumar. 2017. A Review of Multivariate Distributions for Count Data Derived from the Poisson Distribution. Wiley interdisciplinary reviews. Computational statistics 9.
#' @references Raqab, M. Z., D. Kundu, and F. A. Al-Awadhi. 2019. Compound zero-truncated Poisson normal distribution and its applications. Communications in Statistics - Theory and Methods:1–21.
#'
#' @keywords internal
"_PACKAGE"

#' @import Rcpp
#' @import stats
#' @importFrom DistributionUtils is.wholenumber
#' @importFrom mixtools normalmixEM

#' @useDynLib ztpln, .registration=TRUE
NULL
