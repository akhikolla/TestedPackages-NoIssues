\name{FastGaSP-package}
\alias{FastGaSP-package}
\alias{FastGaSP}
\docType{package}
\title{
\packageTitle{FastGaSP}
}
\description{
\packageDescription{FastGaSP}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{FastGaSP}
\packageIndices{FastGaSP}
Fast computational algorithms for Gaussian stochastic process with Matern kernels 
by the forward filtering and backward smoothing algorithm.
}
\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}
\references{

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models},  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.



M. Gu, Y. Xu (2017), \emph{Nonseparable Gaussian stochastic process: a unified
view and computational strategy}, arXiv:1711.11501.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian Stochastic Process Emulation}, \emph{Annals of Statistics}, \bold{46}, 3038-3066.

}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%~~ directory ~~
\keyword{Gaussian stochastic process}
\keyword{Kalman filter}
\keyword{profile likelihood}
\keyword{prediction}


\seealso{
\code{\link[FastGaSP]{FastGaSP}}
}
\examples{

library(FastGaSP)

#------------------------------------------------------------------------------
# Example 1 : fast computation algorithm for noisy data
#------------------------------------------------------------------------------

y_R<-function(x){
  sin(2*pi*x)
}

###let's test for 2000 observations
set.seed(1)
num_obs=2000
input=runif(num_obs)
output=y_R(input)+rnorm(num_obs,mean=0,sd=0.1)

##constucting the fgasp.model
fgasp.model=fgasp(input, output)

##range and noise-variance ratio (nugget) parameters 
param=c( log(1),log(.02))
## the log lik
log_lik(param,fgasp.model)
##time cost to compute the likelihood
time_cost=system.time(log_lik(param,fgasp.model))
time_cost[1]

##consider a nonparametric regression setting 
##estimate the parameter by maximum likelihood estimation

est_all<-optim(c(log(1),log(.02)),log_lik,object=fgasp.model,method="L-BFGS-B",
              control = list(fnscale=-1))  

##estimated log inverse range parameter and log nugget
est_all$par

##estimate variance 
est.var=Get_log_det_S2(est_all$par,fgasp.model@have_noise,fgasp.model@delta_x,
                          fgasp.model@output,fgasp.model@kernel_type)[[2]]/fgasp.model@num_obs
est.var

###1. Do some interpolation test 
num_test=5000
testing_input=runif(num_test) ##there are the input where you don't have observations
pred.model=predict(param=est_all$par,object=fgasp.model,testing_input=testing_input)

lb=pred.model@mean+qnorm(0.025)*sqrt(pred.model@var)
ub=pred.model@mean+qnorm(0.975)*sqrt(pred.model@var)

## calculate lb for the mean function
pred.model2=predict(param=est_all$par,object=fgasp.model,testing_input=testing_input,var_data=FALSE)
lb_mean_funct=pred.model2@mean+qnorm(0.025)*sqrt(pred.model2@var)
ub_mean_funct=pred.model2@mean+qnorm(0.975)*sqrt(pred.model2@var)

## plot the prediction
min_val=min(lb,output)
max_val=max(ub,output)

plot(pred.model@testing_input,pred.model@mean,type='l',col='blue',
     ylim=c(min_val,max_val),
     xlab='x',ylab='y')
polygon(c(pred.model@testing_input,rev(pred.model@testing_input)),
        c(lb,rev(ub)),col = "grey80", border = FALSE)
lines(pred.model@testing_input,pred.model@mean,type='l',col='blue')
lines(pred.model@testing_input,y_R(pred.model@testing_input),type='l',col='black')
lines(pred.model2@testing_input,lb_mean_funct,col='blue',lty=2)
lines(pred.model2@testing_input,ub_mean_funct,col='blue',lty=2)
lines(input,output,type='p',pch=16,col='black',cex=0.4) #one can plot data

legend("bottomleft", legend=c("predictive mean","95\% predictive interval","truth"),
       col=c("blue","blue","black"), lty=c(1,2,1), cex=.8)

#--------------------------------------------------------------
# Example 2: example that one does not have a noise in the data
#--------------------------------------------------------------
## Here is a function in the Sobolev Space with order 3
y_R<-function(x){
  j_seq=seq(1,200,1)
  record_y_R=0
  for(i_j in 1:200){
    record_y_R=record_y_R+2*j_seq[i_j]^{-2*3}*sin(j_seq[i_j])*cos(pi*(j_seq[i_j]-0.5)*x)

  }
  record_y_R
}


##generate some data without noise
num_obs=50
input=seq(0,1,1/(num_obs-1))

output=y_R(input)


##constucting the fgasp.model
fgasp.model=fgasp(input, output,have_noise=FALSE)

##range and noise-variance ratio (nugget) parameters 
param=c( log(1))
## the log lik
log_lik(param,fgasp.model)



#if one does not have noise one may need to give a lower bound or use a penalty 
#(e.g. induced by a prior) to make the estimation more robust
est_all<-optimize(log_lik,interval=c(0,10),maximum=TRUE,fgasp.model)
  
##Do some interpolation test for comparison
num_test=1000
testing_input=runif(num_test) ##there are the input where you don't have observations

pred.model=predict(param=est_all$maximum,object=fgasp.model,testing_input=testing_input)



#This is the 95 posterior credible interval for the outcomes which contain the estimated 
#variance of the noise
#sometimes there are numerical instability is one does not have noise or error
lb=pred.model@mean+qnorm(0.025)*sqrt(abs(pred.model@var))
ub=pred.model@mean+qnorm(0.975)*sqrt(abs(pred.model@var))

## plot the prediction
min_val=min(lb,output)
max_val=max(ub,output)

plot(pred.model@testing_input,pred.model@mean,type='l',col='blue',
     ylim=c(min_val,max_val),
     xlab='x',ylab='y')
polygon( c(pred.model@testing_input,rev(pred.model@testing_input)),
         c(lb,rev(ub)),col = "grey80", border = FALSE)
lines(pred.model@testing_input,pred.model@mean,type='l',col='blue')
lines(pred.model@testing_input,y_R(pred.model@testing_input),type='l',col='black')
lines(input,output,type='p',pch=16,col='black')
legend("bottomleft", legend=c("predictive mean","95\% predictive interval","truth"),
       col=c("blue","blue","black"), lty=c(1,2,1), cex=.8)


##mean square error for all inputs
mean((pred.model@mean- y_R(pred.model@testing_input))^2)

}
