#' @title  Summarizing Initial Sample Size Estimates
#' @description \code{summary} method for class "\code{ssest}".
#'
#' @param object an object of class "\code{ssest}".
#' @param ... Arguments to be passed to or from other methods.
#'
#' @details
#' \code{summary.ssest} gives back initial sample size estimates required. Furthermore, inputs are displayed for double checking.
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
#' summary(estimate)
#' @export

summary.ssest<-function(object, ...){

  if(object$model=="NB-INAR(1)"){
    cat("Initial Sample Size Calculation")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("alpha level:          ", object$alpha))
    cat("\n")
    cat(paste("testing power:        ", round(object$power,2)))
    cat("\n")
    cat(paste("rate ratio:           ", object$delta))
    cat("\n")
    cat(paste("rate control group    ", object$muC))
    cat("\n")
    cat(paste("dispersion parameter: ", object$size))
    cat("\n")
    cat(paste("correlation parameter:", object$rho))
    cat("\n")
    cat(paste("time points:          ", object$tp))
    cat("\n")
    cat(paste("allocation factor:    ", object$k))
    cat("\n\n")
    cat("Sample Size")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("control group:  ", round(object$n[1],2)))
    cat("\n")
    cat(paste("treatment group:"), round(object$n[2],2))
  }

  if(object$model=="GF"){
    cat("Initial Sample Size Calculation")
    cat("\n")
    cat("---------------------------------")
    cat("\n")
    cat(paste("alpha level:          ", object$alpha))
    cat("\n")
    cat(paste("testing power:        ", round(object$power,2)))
    cat("\n")
    cat(paste("trend type:           ", object$trend))
    cat("\n")
    cat(paste(c("trend parameters:     ", object$lambda)))
    cat("\n")
    cat(paste("dispersion parameter: ", object$size))
    cat("\n")
    cat(paste("correlation parameter:", object$rho))
    cat("\n")
    cat(paste("time points:          ", object$tp))
    cat("\n")
    cat(paste("allocation factor:    ", object$k))
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
    cat("Initial Sample Size Calculation")
    cat("\n")
    cat("Model: One Subgroup within a Full Population")
    cat("\n")
    cat("---------------------------------------------")
    cat("\n")
    cat(paste("alpha level:            ", object$alpha))
    cat("\n")
    cat(paste("1-(testing power):      ", object$beta))
    cat("\n")
    cat(paste("prevalence of subgroup: ", object$tau))
    cat("\n")
    cat(paste("effect in subgroup:     ", object$delta[2]))
    cat("\n")
    cat(paste("effect outside subgroup:", object$delta[1]))
    cat("\n")
    cat(paste("sd in subgroup:         ", object$sigma[2]))
    cat("\n")
    cat(paste("sd outside subgroup:    ", object$sigma[1]))
    cat("\n")
    cat(paste("power precision:        ", object$eps))
    cat("\n")
    cat(paste("approximation method:   ", object$approx))
    cat("\n")
    cat(paste("allocation factor:      ", object$k))
    cat("\n\n")
    cat("Sample Size")
    cat("\n")
    cat("---------------------------------------------")
    cat("\n")
    cat(paste("control group:  ", round(object$n[1],2)))
    cat("\n")
    cat(paste("treatment group:"), round(object$n[2],2))
  }
}
