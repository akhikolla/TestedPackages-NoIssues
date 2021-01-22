##' Create monotonic polynomial graded response (GR-MP) model
##'
##' The GR-MP model replaces the linear predictor of the graded
##' response model (Samejima, 1969, 1972) with a monotonic polynomial (Falk, conditionally accepted).
##'
##' Given its relationship to the graded response model, the GR-MP
##' is constructed in an analogous way:
##'
##' \deqn{\mathrm P(\mathrm{pick}=0|\lambda,\alpha,\tau,\theta) = 1- \frac{1}{1+\exp(-(\xi_1 + m(\theta;\lambda,\mathbf{\alpha},\mathbf{\tau})))}
##' }{P(pick=0|a,c,theta) = 1-1/(1+exp(-(xi_1 + m(theta;lambda,alpha,tau)))}
##' \deqn{\mathrm P(\mathrm{pick}=k|\lambda,\alpha,\tau,\theta) = \frac{1}{1+\exp(-(\xi_k + m(\theta;\lambda,\mathbf{\alpha},\mathbf{\tau})))} - \frac{1}{1+\exp(-(\xi_{k+1} + m(\theta,\lambda,\mathbf{\alpha},\mathbf{\tau}))}
##' }{P(pick=k|a,c,theta) = 1/(1+exp(-(xi_k + m(theta;lambda,alpha,tau))) - 1/(1+exp(-(xi_(k+1) + m(theta;lambda,alpha,tau))))}
##' \deqn{\mathrm P(\mathrm{pick}=K|\lambda,\alpha,\tau,\theta) = \frac{1}{1+\exp(-(\xi_K + m(\theta;\lambda,\mathbf{\alpha},\mathbf{\tau}))}
##' }{P(pick=K|\lambda,\alpha,\tau,theta) = 1/(1+exp(-(xi_K + m(theta;lambda,alpha,tau))))}
##'
##' The order of the polynomial is always odd and is controlled by
##' the user specified non-negative integer, q. The model contains
##' 1+(outcomtes-1)+2*q parameters and are used as input to the \code{\link{rpf.prob}}
##' or \code{\link{rpf.dTheta}} functions in the following order:
##' \eqn{\lambda}{lambda} - slope of the item model when q=0,
##' \eqn{\xi}{xi} - a (outcomes-1)-length vector of intercept parameters,
##' \eqn{\alpha}{alpha} and \eqn{\tau}{tau} - two parameters that control bends in
##' the polynomial. These latter parameters are repeated in the same order for
##' models with q>0.  For example, a q=2 polynomial with 3 categories will have an item
##' parameter vector of: \eqn{\lambda, \xi_1, \xi_2, \alpha_1, \tau_1, \alpha_2, \tau_2}{
##' lambda, xi1, xi2, alpha1, tau1, alpha2, tau2}.
##'
##' As with other monotonic polynomial-based item models
##' (e.g., \code{\link{rpf.lmp}}), the polynomial looks like the
##' following:
##'
##' \deqn{m(\theta;\lambda,\alpha,\tau) = b_1\theta + b_2\theta^2 + \dots + b_{2q+1}\theta^{2q+1}
##' }{m(theta) = b_1*theta + b_2*theta^2 + \dots + b_(2q+1)*theta^{2q+1}}
##'
##' However, the coefficients, b, are not directly estimated, but are a function of the
##' item parameters, and the parameterization of the GR-MP is different than
##' that currently appearing for the logistic function of a monotonic
##' polynomial (LMP; \code{\link{rpf.lmp}}) and monotonic polynomial generalized partial credit
##' (GPC-MP; \code{\link{rpf.gpcmp}}) models. In particular, the polynomial is
##' parameterized such that boundary descrimination functions for the GR-MP will
##' be all monotonically increasing or decreasing for any given item. This allows
##' the possibility of items that load either negatively or positively on the latent
##' trait, as is common with reverse-worded items in non-cognitive tests (e.g., personality).
##'
##' The derivative \eqn{m'(\theta;\lambda,\alpha,\tau)}{m'(theta;lambda,alpha,tau)} is
##' parameterized in the following way:
##'
##'\deqn{ m'(\theta;\lambda,\alpha,\tau) = \left\{\begin{array}{ll}\lambda \prod_{u=1}^q(1-2\alpha_{u}\theta + (\alpha_{u}^2 + \exp(\tau_{u}))\theta^2) &  \mbox{if } q > 0 \\
##'  \lambda & \mbox{if } q = 0\end{array} \right.}{m'(theta) = m'(theta;lambda,alpha,tau) = lambda \prod_{u=1}^q (1-2*alpha_u*theta + (alpha_u^2 + exp(tau_u))*theta^2) (if q > 0) \\
##'  lambda (if q = 0)}
##'
##' Note that the only difference between the GR-MP and these other models
##' is that \eqn{\lambda}{lambda} is not re-parameterized and may take on
##' negative values. When \eqn{\lambda}{lambda} is negative, it is analogous
##' to having a negative loading or a monotonically decreasing function.
##'
##' @param outcomes The number of possible response categories. When equal to 2, the model reduces
##' to the logistic function of a monotonic polynomial (LMP).
##' @param q a non-negative integer that controls the order of the
##' polynomial (2q+1) with a default of q=0 (1st order polynomial = graded response model).
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{FALSE}. The multidimensional version is not yet
##' available.
##' @return an item model
##' @references Falk, C. F. (conditionally accepted). The monotonic polynomial
##' graded response model: Implementation and a comparative study. \emph{Applied Psychological
##' Measurement}.
##' @references Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. \emph{Psychometric Monographs}, 17.
##' @references Samejima, F. (1972). A general model of free-response data. \emph{Psychometric Monographs}, 18.
##'
##' @family response model
##'
##' @examples
##' spec <- rpf.grmp(5,2) # 5-category, 3rd order polynomial
##' theta<-seq(-3,3,.1)
##' p<-rpf.prob(spec, c(2.77,2,1,0,-1,.89,-8.7,-.74,-8.99),theta)

rpf.grmp <- function(outcomes=2, q=0, multidimensional=FALSE) {
  if(!(q%%1==0)){
    stop("q must be an integer >= 0")
  }
  if(multidimensional){
      stop("Multidimensional grmp model is not yet supported")
  }
  m <- NULL
  id <- -1
  id <- rpf.id_of("grmp")
  m <- new("rpf.1dim.grmp",
           outcomes=outcomes,
           factors=1)
  m@spec <- c(id, m@outcomes, m@factors, q)
  m
}

setMethod("rpf.rparam", signature(m="rpf.1dim.grmp"),
          function(m, version) {
            n <- 1
            q<-m$spec[4] ## ok to hardcode this index?
            ret<-c(lambda=rlnorm(n, 0, .5)) # random overall slope
            # randomly generate xi
            xi <- rnorm(m@outcomes-1)
            xi <- xi[order(-xi)]
            ret<-c(ret,xi=xi)
            if(q>0){
                for(i in 1:q){
                    ret<-c(ret,runif(n,-1,1),log(runif(n,.0001,1)))
                    names(ret)[(m@outcomes+1+(i-1)*2):(m@outcomes+(i*2))]<-c(paste("alpha",i,sep=""),paste("tau",i,sep=""))
                }
            }
            ret
        })
