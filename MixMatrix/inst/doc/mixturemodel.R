## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(MixMatrix)

## ----demo----------------------------------------------------------------
library(MixMatrix)
 set.seed(20180221)
 A <- rmatrixt(30,mean=matrix(0,nrow=3,ncol=4), df = 10) # 3x4 matrices with mean 0
 B <- rmatrixt(30,mean=matrix(1,nrow=3,ncol=4), df = 10) # 3x4 matrices with mean 2
 C <- array(c(A,B), dim=c(3,4,60)) # combine into one array
 prior <- c(.5,.5) # equal probability prior
 # create an intialization object, starts at the true parameters
 init = list(centers = array(c(rep(0,12),rep(1,12)), dim = c(3,4,2)),
              U = array(c(diag(3), diag(3)), dim = c(3,3,2)),
              V = array(c(diag(4), diag(4)), dim = c(4,4,2))
             )
 # fit model
 res<-matrixmixture(C, init = init, prior = prior, nu = 10,
                    model = "t", tolerance = 1e-2)
 print(res$centers) # the final centers
 print(res$pi) # the final mixing proportion
 logLik(res)
 AIC(logLik(res))
 plot(res) # the log likelihood by iteration


## ----initializer---------------------------------------------------------

init_matrixmixture(C, prior = c(.5,.5), centermethod = 'kmeans')

init_matrixmixture(C, K = 2, centermethod = 'random')


## ----final---------------------------------------------------------------
sessionInfo()


## ----getlabels, echo = FALSE---------------------------------------------
labs = knitr::all_labels()
labs = labs[!labs %in% c("setup", "toc", "getlabels", "allcode")]

## ----allcode, ref.label = labs, eval = FALSE-----------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>"
#  )
#  library(MixMatrix)
#   set.seed(20180221)
#   A <- rmatrixt(30,mean=matrix(0,nrow=3,ncol=4), df = 10) # 3x4 matrices with mean 0
#   B <- rmatrixt(30,mean=matrix(1,nrow=3,ncol=4), df = 10) # 3x4 matrices with mean 2
#   C <- array(c(A,B), dim=c(3,4,60)) # combine into one array
#   prior <- c(.5,.5) # equal probability prior
#   # create an intialization object, starts at the true parameters
#   init = list(centers = array(c(rep(0,12),rep(1,12)), dim = c(3,4,2)),
#                U = array(c(diag(3), diag(3)), dim = c(3,3,2)),
#                V = array(c(diag(4), diag(4)), dim = c(4,4,2))
#               )
#   # fit model
#   res<-matrixmixture(C, init = init, prior = prior, nu = 10,
#                      model = "t", tolerance = 1e-2)
#   print(res$centers) # the final centers
#   print(res$pi) # the final mixing proportion
#   logLik(res)
#   AIC(logLik(res))
#   plot(res) # the log likelihood by iteration
#  
#  
#  init_matrixmixture(C, prior = c(.5,.5), centermethod = 'kmeans')
#  
#  init_matrixmixture(C, K = 2, centermethod = 'random')
#  
#  sessionInfo()
#  

