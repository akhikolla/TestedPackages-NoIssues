## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(MixMatrix)

## ----generatedata--------------------------------------------------------
set.seed(20180222)
library('MixMatrix')
A <- rmatrixnorm(30, mean = matrix(0, nrow=2, ncol=3)) # creating the three groups
B <- rmatrixnorm(30, mean = matrix(c(1,0), nrow=2, ncol=3))
C <- rmatrixnorm(30, mean = matrix(c(0,1), nrow=2, ncol=3))
ABC <- array(c(A,B,C), dim = c(2,3,90)) # combining into on array
groups <- factor(c(rep("A",30),rep("B",30),rep("C",30))) # labels
prior = c(30,30,30)/90 # equal prior
matlda <- matrixlda(x = ABC,grouping = groups, prior = prior) # perform LDA
matqda <- matrixqda(x = ABC,grouping = groups, prior = prior) # perform QDA

## ----predict-------------------------------------------------------------
ABC[,,c(1,31,61)] # true class memberships: A, B, C
#predict the membership of the first observation of each group
predict(matlda, ABC[,,c(1,31,61)]) 
#predict the membership of the first observation of each group
predict(matqda,  ABC[,,c(1,31,61)])


## ----objectstructure-----------------------------------------------------
matlda

matqda


## ----sessioninfo---------------------------------------------------------
sessionInfo()

## ----getlabels, echo = FALSE---------------------------------------------
labs = knitr::all_labels()
labs = labs[!labs %in% c("setup", "toc", "getlabels", "allcode")]

## ----allcode, ref.label = labs, eval = FALSE-----------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>"
#  )
#  set.seed(20180222)
#  library('MixMatrix')
#  A <- rmatrixnorm(30, mean = matrix(0, nrow=2, ncol=3)) # creating the three groups
#  B <- rmatrixnorm(30, mean = matrix(c(1,0), nrow=2, ncol=3))
#  C <- rmatrixnorm(30, mean = matrix(c(0,1), nrow=2, ncol=3))
#  ABC <- array(c(A,B,C), dim = c(2,3,90)) # combining into on array
#  groups <- factor(c(rep("A",30),rep("B",30),rep("C",30))) # labels
#  prior = c(30,30,30)/90 # equal prior
#  matlda <- matrixlda(x = ABC,grouping = groups, prior = prior) # perform LDA
#  matqda <- matrixqda(x = ABC,grouping = groups, prior = prior) # perform QDA
#  ABC[,,c(1,31,61)] # true class memberships: A, B, C
#  #predict the membership of the first observation of each group
#  predict(matlda, ABC[,,c(1,31,61)])
#  #predict the membership of the first observation of each group
#  predict(matqda,  ABC[,,c(1,31,61)])
#  
#  matlda
#  
#  matqda
#  
#  sessionInfo()

