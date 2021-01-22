#' @title  Summarizing Blinded Sample Size Reestimation
#' @description \code{summary} method for class "\code{bssrest}".
#'
#' @param object an object of class "\code{bssrest}".
#' @param ... Arguments to be passed to or from other methods.
#'
#' @details
#' \code{summary.bssrest} gives back blinded sample size estimates. Furthermore, inputs are displayed for double checking.
#'
#' @seealso \code{\link{n.nb.inar1}} for initial sample size estimates within the NB-INAR(1) model.
#'
#' @examples
#' #Calculate required sample size to find significant difference with
#' #80% probability when testing the Nullhypothesis H_0: mu_T/mu_C >= 1
#' #assuming the true effect delta is 0.8 and rate, size and correlation
#' #parameter in the control group are 2, 1 and 0.5, respectively.
#'
#' estimate<-n.nb.inar1(alpha=0.025, power=0.8, delta=0.8, muC=2, size=1, rho=0.5, tp=7, k=1)
#'
#' #Simulate data
#' set.seed(8)
#' placebo<-rnbinom.inar1(n=50, size=1, mu=2, rho=0.5, tp=7)
#' treatment<-rnbinom.inar1(n=50, size=1, mu=1.6, rho=0.5, tp=7)
#'
#' #Blinded sample size reestimation
#' estimate<-bssr.nb.inar1(alpha=0.025, power=0.8, delta=0.8, x=rbind(placebo, treatment),
#'   n=c(50,50), k=1)
#' summary(estimate)
#' @export

summary.bssrest<-function(object, ...){

  if(object$model=="NB-INAR(1)"){
    cat("Blinded Sample Size Reestimation")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("alpha level:               ", object$alpha))
    cat("\n")
    cat(paste("testing power:             ", object$beta))
    cat("\n")
    cat(paste("rate ratio:                ", object$delta))
    cat("\n")
    cat(paste("est. overall rate          ", round(object$mu, 2)))
    cat("\n")
    cat(paste("est. dispersion parameter: ", round(object$size, 2)))
    cat("\n")
    cat(paste("est. correlation parameter:", round(object$rho, 2)))
    cat("\n")
    cat(paste("time points:               ", object$tp))
    cat("\n")
    cat(paste("desired allocation factor: ", object$k))
    cat("\n\n")
    cat("Reestimated Sample Size")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("control group:  ", round(object$n[1],2)))
    cat("\n")
    cat(paste("treatment group:"), round(object$n[2],2))
  }

  if(object$model=="GF"){
    cat("Blinded Sample Size Reestimation")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("alpha level:               ", object$alpha))
    cat("\n")
    cat(paste("testing power:             ", round(object$power,2)))
    cat("\n")
    cat(paste("trend type:                ", object$trend))
    cat("\n")
    cat(paste("effect size:               ", object$delta))
    cat("\n")
    cat(paste(c("est. trend parameters:     ", round(object$lambda,2))))
    cat("\n")
    cat(paste("est. dispersion parameter: ", round(object$size,2)))
    cat("\n")
    cat(paste("est. correlation parameter:", round(object$rho,2)))
    cat("\n")
    cat(paste("time points:               ", object$tp))
    cat("\n")
    cat(paste("allocation factor:         ", object$k))
    cat("\n\n")
    cat("Sample Size")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("control group:  ", round(object$n[1],2)))
    cat("\n")
    cat(paste("treatment group:"), round(object$n[2],2))
  }

  if(object$model=="normal1subgroup"){
    cat("Blinded Sample Size Reestimation")
    cat("\n")
    cat("Model: One Subgroup within a Full Population")
    cat("\n")
    cat("-----------------------------------------------")
    cat("\n")
    cat(paste("alpha level:               ", object$alpha))
    cat("\n")
    cat(paste("1-(testing power):         ", object$beta))
    cat("\n")
    cat(paste("effect in subgroup:        ", object$delta[2]))
    cat("\n")
    cat(paste("effect outside subgroup:   ", object$delta[1]))
    cat("\n")
    cat(paste("est. prevalence subgroup:  ", round(object$tau.est,3)))
    cat("\n")
    cat(paste("est. sd in subgroup:       ", round(object$sigma.est[2],3)))
    cat("\n")
    cat(paste("est. sd outside subgroup:  ", round(object$sigma.est[1],3)))
    cat("\n")
    cat(paste("power precision:           ", object$eps))
    cat("\n")
    cat(paste("approximation method:      ", object$approx))
    cat("\n")
    cat(paste("allocation factor:         ", object$k))
    cat("\n\n")
    cat("Reestimated Sample Size")
    cat("\n")
    cat("-----------------------------------------------")
    cat("\n")
    cat(paste("control group:  ", round(object$n[1],2)))
    cat("\n")
    cat(paste("treatment group:"), round(object$n[2],2))
  }
}
