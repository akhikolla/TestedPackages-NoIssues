##################################################
### Percentile-Calibrated ('perc-cal') Confidence Interval Functions

library(RcppEigen)
library(Rcpp)
library(RcppArmadillo)

##' Calculates Percentile-Calibrated Linear Regression Confidence Intervals
##'
##' This is the main function of the package.  It takes as inputs the 
##' predictor/response matrix appended together, which can be either
##' a data frame or a matrix, along with the desired coverage and other
##' settings, and outputs marginal confidence intervals for each of the 
##' predictors, including the intercept.
##' @name perccal_interval
##' @param Xy [n by (p+1)] matrix: X in columns 1 to p, y in column p+1. 
##' X is the design matrix, and is assumed to not include a vector of one's.
##' @param alpha Target coverage desired.
##' @param G Number of grid points to evaluate calibrated percentile method 
##' on each side over.
##' @param B Number of 1st stage bootstrap samples.
##' @param B2 Number of 2nd stage double bootstrap samples.
##' @return Return a (p+1)x2 matrix containing confidence intervals for all 
##' regression coefficients, estimated via the perc-cal method.
##' @examples 
##' set.seed(1234)
##' n = 32
##' B = 500
##' B2 = 500
##' G=20
##' x1=rnorm(n)
##' x2=rnorm(n)
##' eps=rnorm(n)
##' y = x1 + 2*x2 + eps
##' Xy = cbind(x1,x2,y)
##' alpha = .025	
##' perccal_interval(Xy, alpha, G, B, B2)
perccal_interval  = function(Xy, alpha, G = 20, B = 999, B2 = 999){
	# Inputs: 
	# Xy: [n by (p+1)] matrix: X in columns 1 to p, y in column p+1 
				#(please do not include the intercept)
	# alpha: Target coverage desired.
	# G    : # of grid points to evaluate calibrated percentile method on each side.
	# B    : # of 1st stage bootstrap samples
	# B2   : # of 2nd stage bootstrap samples
	# Output: 
	# Confidence interval for all regression coefficients, estimated via perc-cal method.
	
	# Initialize all the matrices and vectors to be used below.
	p = ncol(Xy)-1
	theta.hat.boot = matrix(NA,nrow=p+1,ncol=B)
	theta.hat.dboot = matrix(NA,nrow=p+1,ncol=B2)
	theta.qtl.lgrid.lo = array(NA, dim=c(p+1,B,G))
	theta.qtl.lgrid.hi = array(NA, dim=c(p+1,B,G))
	perc.cal.ints = matrix(NA,nrow=p+1,ncol=2)
	
	l.grid.lo = seq(.0001, alpha*1.4, length.out=G) 
	l.grid.hi = 1-l.grid.lo
		
	# Start main code
	X=as.matrix(Xy[,1:p])
	Xm = cbind(1,X)
	y=Xy[,p+1]
	regr_output = coefficients(summary(lm(y~X)))
	prednames = rownames(regr_output)
	theta.hat = as.numeric(regr_output[,"Estimate"])
	
	boot_output = Cdboot_multi(as.matrix(Xy), l.grid.lo, l.grid.hi, B, B2, G)
	
	for(k in 1:(p+1)){
		theta.hat.boot	= boot_output$theta_hat_boot[,k]
		theta.qtl.lgrid.lo = boot_output$theta_qtl_lgrid_lo[k,1][[1]][,,1]
		theta.qtl.lgrid.hi = boot_output$theta_qtl_lgrid_hi[k,1][[1]][,,1]
		
		# perc-cal method
		log.lgrid.lo = theta.qtl.lgrid.lo < theta.hat[k]  
		log.lgrid.hi = theta.qtl.lgrid.hi > theta.hat[k]        
		percs.lo = 1-colSums(log.lgrid.lo)/B
		percs.hi = 1-colSums(log.lgrid.hi)/B
		lo.log = (percs.lo[1:(G-1)]<alpha)*(percs.lo[2:G]>=alpha)
		hi.log = (percs.hi[1:(G-1)]<alpha)*(percs.hi[2:G]>=alpha)
		if(sum(lo.log)>0){l.alpha.lo = l.grid.lo[lo.log*(1:(G-1))]}
		if(sum(lo.log)<=0){l.alpha.lo = .0001}
		if(sum(hi.log)>0){l.alpha.hi = l.grid.hi[hi.log*(1:(G-1))]}
		if(sum(hi.log)<=0){l.alpha.hi = 1-.0001}
		perc.cal.ints[k,] = quantile(theta.hat.boot,c(l.alpha.lo, l.alpha.hi))
		}
	
	colnames(perc.cal.ints) = c(paste("lower ",round(alpha,3),sep=""),paste("upper ",1-round(alpha,3),sep=""))
	rownames(perc.cal.ints) = prednames
	perc.cal.ints
	}


	


##' Fast computation of quantiles
##'
##' Helper function which takes as input a vector and obtains quantiles for it.
##' Number of quantiles may be greater than one.
##' @name Cquantile
##' @param xx Numeric vector we are obtaining quantiles for.
##' @param p Numeric vector of quantiles.
##' @return Numeric vector containing quantiles, possibly greater than one.
Cquantile <- function(xx, p) {
    .Call('perccal_Cquantile', PACKAGE = 'perccal', xx, p)
}


##' Sample from [1,2,...,N] with replacement nsamp times in Rcpp.
##'
##' Helper function which samples from [1,2,...,N] with replacement 
##' nsamp times in Rcpp.
##' @name sample_rcpp
##' @param N Largest integer to sample from.
##' @param nsamp number of samples from [1,2,...,N] with replacement to obtain.
##' @return samps nsamp-length vector of samples from [1,2,...,N] with replacement to obtain.
sample_rcpp <- function(N, nsamp) {
    .Call('perccal_sample_rcpp', PACKAGE = 'perccal', N, nsamp)
}
	
	

##' Fast computation of internal double bootstrap calculations
##'
##' This is the workhorse function of the package, speeding up computations
##' within double bootstrap routine.
##' @name Cdboot_multi
##' @param xxyy (n by p+1) matrix for X (design matrix) and response vector y.
##' @param lgridlo Lower quantile values of double bootstrap distribution to 
##' obtain.
##' @param lgridhi Upper quantile values of double bootstrap distribution to 
##' obtain.
##' @param B Number of 1st stage bootstrap samples.
##' @param B2 Number of 2nd stage double bootstrap samples.
##' @param G Calculate quantile-based empirical coverage at this many grid 
##' points
##' @return theta_hat_boot first-level bootstrap estimates of all slope coefficients
##' @return theta_qtl_lgrid_lo (p+1 by B by G by 1) matrix for lower quantiles at 
##' all grid points for all predictors over all bootstrap samples.
##' @return theta_qtl_lgrid_hi (p+1 by B by G by 1) matrix for upper quantiles at 
##' all grid points for all predictors over all bootstrap samples.
Cdboot_multi <- function(xxyy, lgridlo, lgridhi, B, B2, G) {
    .Call('perccal_Cdboot_multi', PACKAGE = 'perccal', xxyy, lgridlo, lgridhi, B, B2, G)
}

	




	
	
	
	
	

