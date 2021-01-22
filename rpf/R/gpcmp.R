##' Create monotonic polynomial generalized partial credit (GPC-MP) model
##'
##' This model is a polytomous model proposed by Falk & Cai (2016) and is based
##' on the generalized partial credit model (Muraki, 1992).
##'
##' The GPC-MP replaces the linear predictor part of the
##' generalized partial credit model with a monotonic polynomial,
##' \eqn{m(\theta;\omega,\xi,\mathbf{\alpha},\mathbf{\tau})}{m(theta; omega, alpha, tau)}.
##' The response function for category k is:
##'
##' \deqn{\mathrm P(\mathrm{pick}=k|\omega,\xi,\alpha,\tau,\theta)
##' = \frac{\exp(\sum_{v=0}^k (\xi_k + m(\theta;\omega,\xi,\mathbf{\alpha},\mathbf{\tau})))}{\sum_{u=0}^{K-1}\exp(\sum_{v=0}^u (\xi_u + m(\theta;\omega,\xi,\mathbf{\alpha},\mathbf{\tau})))}
##' }{\frac{exp(sum_{v=0}^k (xi_k + m(theta;omega,xi,alpha,tau)))}{sum_{u=0}^{K-1}exp(sum_{v=0}^u (xi_u + m(theta;omega,xi,alpha,tau)))}}
##'
##' where \eqn{\mathbf{\alpha}}{alpha} and \eqn{\mathbf{\tau}}{tau} are vectors
##' of length q. The GPC-MP uses the same parameterization for the polynomial
##' as described for the logistic function of a monotonic polynomial (LMP).
##' See also (\code{\link{rpf.lmp}}).
##'
##' The order of the polynomial is always odd and is controlled by
##' the user specified non-negative integer, q. The model contains
##' 1+(outcomtes-1)+2*q parameters and are used as input to the \code{\link{rpf.prob}}
##' function in the following order:
##' \eqn{\omega}{omega} - natural log of the slope of the item model when q=0,
##' \eqn{\xi}{xi} - a (outcomes-1)-length vector of intercept parameters,
##' \eqn{\alpha}{alpha} and \eqn{\tau}{tau} - two parameters that control bends in
##' the polynomial. These latter parameters are repeated in the same order for
##' models with q>0. For example, a q=2 polynomial with 3 categories will have an item
##' parameter vector of: \eqn{\omega, \xi_1, \xi_2, \alpha_1, \tau_1, \alpha_2, \tau_2}{
##' omega, xi1, xi2, alpha1, tau1, alpha2, tau2}.
##'
##' Note that the GPC-MP reduces to the LMP when the number of
##' categories is 2, and the GPC-MP reduces to the generalized partial
##' credit model when the order of the polynomial is 1 (i.e., q=0).
##'
##' @param outcomes The number of possible response categories.
##' @param q a non-negative integer that controls the order of the
##' polynomial (2q+1) with a default of q=0 (1st order polynomial = generalized partial credit model).
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{FALSE}. The multidimensional version is not yet
##' available.
##' @return an item model
##' @references Falk, C. F., & Cai, L. (2016). Maximum marginal likelihood
##' estimation of a monotonic polynomial generalized partial credit model with
##' applications to multiple group analysis. \emph{Psychometrika, 81}, 434-460.
##' \url{http://dx.doi.org/10.1007/s11336-014-9428-7}
##' @references Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. \emph{Applied Psychological Measurement, 16,} 159–176.
##'
##' @family response model
##'
##' @examples
##' spec <- rpf.gpcmp(5,2) # 5-category, 3rd order polynomial
##' theta<-seq(-3,3,.1)
##' p<-rpf.prob(spec, c(1.02,3.48,2.5,-.25,-1.64,.89,-8.7,-.74,-8.99),theta)
##'

rpf.gpcmp <- function(outcomes=2, q=0, multidimensional=FALSE) {
  if(!(q%%1==0)){
    stop("q must be an integer >= 0")
  }
  if(multidimensional){
      stop("Multidimensional gpcmp model is not yet supported")
  }
  m <- NULL
  id <- -1
  id <- rpf.id_of("gpcmp")
  m <- new("rpf.1dim.gpcmp",
           outcomes=outcomes,
           factors=1)
  m@spec <- c(id, m@outcomes, m@factors, q)
  m
}

setMethod("rpf.rparam", signature(m="rpf.1dim.gpcmp"),
          function(m, version) {
            n <- 1
            q<-m$spec[4] ## ok to hardcode this index?
            ret<-c(omega=rnorm(n, 0, .5)) # random overall slope
            # randomly generate xi
            xi <- rnorm(m@outcomes-1)
            xi <- xi[order(xi)]
            ret<-c(ret,xi=xi)
            if(q>0){
                for(i in 1:q){
                    ret<-c(ret,runif(n,-1,1),log(runif(n,.0001,1)))
                    names(ret)[(m@outcomes+1+(i-1)*2):(m@outcomes+(i*2))]<-c(paste("alpha",i,sep=""),paste("tau",i,sep=""))
                }
            }
            ret
        })
