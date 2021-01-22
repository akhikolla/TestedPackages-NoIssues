## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(MixMatrix)

## ----mleone--------------------------------------------------------------
set.seed(20190622)
sigma = (1/7) * rWishart(1, 7, 1*diag(3:1))[,,1]
A = rmatrixt(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
   V = sigma, df = 7)
results=MLmatrixt(A, df = 7)
print(results)

## ----dontrun, eval=FALSE, echo = FALSE-----------------------------------
#  ### Here is the long simulation
#  library(ggplot2)
#  
#  set.seed(20181102)
#  
#  df = c(5, 10, 20)
#  df5 <- rep(0,200)
#  df10 <- rep(0,200)
#  df100 <- rep(0,200)
#  df550 <-  rep(0,200)
#  df1050 <-  rep(0,200)
#  df2050 <- rep(0,200)
#  df5100 <-  rep(0,200)
#  df10100 <-  rep(0,200)
#  df20100 <- rep(0,200)
#  
#  meanmat = matrix(0,5,3)
#  U = diag(5)
#  V = diag(3)
#  
#  for(i in 1:200){
#  	df5[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 5, n = 35, U =U, V =V), fixed = FALSE)$nu
#  	df10[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 10, n = 35, U =U, V =V), fixed = FALSE)$nu
#  	df100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 20, n = 35, U =U, V =V), fixed = FALSE)$nu
#  	df550[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 5, n = 50, U =U, V =V), fixed = FALSE)$nu
#  	df1050[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 10, n = 50, U =U, V =V), fixed = FALSE)$nu
#  	df2050[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 20, n = 50, U =U, V =V), fixed = FALSE)$nu
#  	df5100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 5, n = 100, U =U, V =V), fixed = FALSE)$nu
#  	df10100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 10, n = 100, U =U, V =V), fixed = FALSE)$nu
#  	df20100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 20, n = 100, U =U, V =V), fixed = FALSE)$nu
#  }
#  
#  truedataframe = data.frame(truedf = factor(c(5,10,20),
#  	                                       label = c('5 df', '10 df', '20 df')),
#  	                       estdf = c(5,10,20))
#  
#  dfdataframe = data.frame(truedf = factor(rep(rep(c(5,10,20), each = 200),3),
#  	                                     label = c('5 df', '10 df', '20 df')),
#                           estdf = c(df5, df10, df100, df550, df1050, df2050, df5100, df10100, df20100),
#                           samplesize = factor(rep(c(35,50,100), each = 600)))
#  library(tidyverse)
#  
#  denseplot <- ggplot(data = subset(dfdataframe, estdf < 200),
#  	                aes(x=estdf, fill=samplesize)) +
#  	geom_density(alpha = .5) +
#      geom_vline(data = truedataframe,
#                 mapping = aes(xintercept = estdf),
#                 size = .5) +
#      theme_bw() +
#      theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),
#  		  strip.text = element_text(size = 8),
#            legend.justification=c(1,0), legend.position=c(.95,.4),
#            legend.background = element_blank(),
#            legend.text =element_text(size = 8), legend.title = element_text(size = 8)) +
#      ggtitle("Density plot of estimated degrees of freedom compared to actual") +
#      xlab(NULL) +
#      ylab(NULL) +
#      scale_fill_manual(values = c("#050505", "#E69F00", "#56B4E9"),
#  		name = "Sample Size") +
#      facet_wrap(factor(truedf)~.,  scales="free") +
#      NULL
#  denseplot
#  
#  
#  knitr::kable(dfdataframe %>% group_by(truedf, samplesize) %>%
#  	summarize(min = min(estdf), max = max(estdf),
#  	median = median(estdf),
#  	mean=mean(estdf),
#  	sd = sd(estdf)))
#  ### Here ends the long simulation
#  

## ----dftenexample, echo = FALSE, message = FALSE, warning = FALSE--------
##### Here is what is really run

set.seed(20190621)
df10 <- rep(0,50)
df1050 <-  rep(0,50)
df10100 <-  rep(0,50)

for(i in 1:50){
   df10[i] = suppressWarnings(MLmatrixt(rmatrixt(mean = matrix(0,5,3),df = 10, n = 25), fixed = FALSE, df = 5, max.iter = 20)$nu)
   df1050[i] = suppressWarnings(MLmatrixt(rmatrixt(mean = matrix(0,5,3),df = 10, n = 50), fixed = FALSE, df = 5, max.iter = 20)$nu)
   df10100[i] = suppressWarnings(MLmatrixt(rmatrixt(mean = matrix(0,5,3),df = 10, n = 100), fixed = FALSE, df = 5, max.iter = 20)$nu)
}


dfdataframe = data.frame(label = c('10 df'),
                         estdf = c(df10, df1050, df10100),
                         samplesize = factor(rep(c(25,50,100), each = 50)))
library(ggplot2)
library(dplyr)
library(magrittr)
denseplot <- ggplot(data = subset(dfdataframe, estdf < 200),aes(x=estdf, fill=samplesize)) +
    geom_density(alpha = .5) +
    geom_vline(mapping = aes(xintercept = 10),
               size = .5) +
    theme_bw() +
    theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), strip.text = element_text(size = 8), 
          legend.justification=c(1,0), legend.position=c(.95,.4),
          legend.background = element_blank(),
          legend.text =element_text(size = 8), legend.title = element_text(size = 8)) +
    ggtitle("Density plot of estimated degrees of freedom compared to actual") +
    xlab(NULL) +
    ylab(NULL) +     
    scale_fill_manual(values = c("#050505", "#E69F00", "#56B4E9"), name = "Sample Size") +
 #   facet_wrap(factor(truedf)~.,  scales="free") +
    NULL
denseplot


knitr::kable(dfdataframe %>% group_by(samplesize) %>% 
	summarize(min = min(estdf), max = max(estdf),
	median = median(estdf),
	mean=mean(estdf),
	sd = sd(estdf)))

#### Here ends what is really run


## ----genlda--------------------------------------------------------------
A <- rmatrixt(30, mean = matrix(0, nrow=2, ncol=3), df = 10)
B <- rmatrixt(30, mean = matrix(c(1,0), nrow=2, ncol=3), df = 10)
C <- rmatrixt(30, mean = matrix(c(0,1), nrow=2, ncol=3), df = 10)
ABC <- array(c(A,B,C), dim = c(2,3,90))
groups <- factor(c(rep("A",30),rep("B",30),rep("C",30)))
prior = c(30,30,30)/90
matlda <- matrixlda(x = ABC,grouping = groups, prior = prior,
                    method = 't', nu = 10, fixed = TRUE)
predict(matlda, newdata = ABC[,,c(1,31,61)])


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
#  set.seed(20190622)
#  sigma = (1/7) * rWishart(1, 7, 1*diag(3:1))[,,1]
#  A = rmatrixt(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#     V = sigma, df = 7)
#  results=MLmatrixt(A, df = 7)
#  print(results)
#  ### Here is the long simulation
#  library(ggplot2)
#  
#  set.seed(20181102)
#  
#  df = c(5, 10, 20)
#  df5 <- rep(0,200)
#  df10 <- rep(0,200)
#  df100 <- rep(0,200)
#  df550 <-  rep(0,200)
#  df1050 <-  rep(0,200)
#  df2050 <- rep(0,200)
#  df5100 <-  rep(0,200)
#  df10100 <-  rep(0,200)
#  df20100 <- rep(0,200)
#  
#  meanmat = matrix(0,5,3)
#  U = diag(5)
#  V = diag(3)
#  
#  for(i in 1:200){
#  	df5[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 5, n = 35, U =U, V =V), fixed = FALSE)$nu
#  	df10[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 10, n = 35, U =U, V =V), fixed = FALSE)$nu
#  	df100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 20, n = 35, U =U, V =V), fixed = FALSE)$nu
#  	df550[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 5, n = 50, U =U, V =V), fixed = FALSE)$nu
#  	df1050[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 10, n = 50, U =U, V =V), fixed = FALSE)$nu
#  	df2050[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 20, n = 50, U =U, V =V), fixed = FALSE)$nu
#  	df5100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 5, n = 100, U =U, V =V), fixed = FALSE)$nu
#  	df10100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 10, n = 100, U =U, V =V), fixed = FALSE)$nu
#  	df20100[i] = MLmatrixt(rmatrixt(mean = meanmat,
#  		df = 20, n = 100, U =U, V =V), fixed = FALSE)$nu
#  }
#  
#  truedataframe = data.frame(truedf = factor(c(5,10,20),
#  	                                       label = c('5 df', '10 df', '20 df')),
#  	                       estdf = c(5,10,20))
#  
#  dfdataframe = data.frame(truedf = factor(rep(rep(c(5,10,20), each = 200),3),
#  	                                     label = c('5 df', '10 df', '20 df')),
#                           estdf = c(df5, df10, df100, df550, df1050, df2050, df5100, df10100, df20100),
#                           samplesize = factor(rep(c(35,50,100), each = 600)))
#  library(tidyverse)
#  
#  denseplot <- ggplot(data = subset(dfdataframe, estdf < 200),
#  	                aes(x=estdf, fill=samplesize)) +
#  	geom_density(alpha = .5) +
#      geom_vline(data = truedataframe,
#                 mapping = aes(xintercept = estdf),
#                 size = .5) +
#      theme_bw() +
#      theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),
#  		  strip.text = element_text(size = 8),
#            legend.justification=c(1,0), legend.position=c(.95,.4),
#            legend.background = element_blank(),
#            legend.text =element_text(size = 8), legend.title = element_text(size = 8)) +
#      ggtitle("Density plot of estimated degrees of freedom compared to actual") +
#      xlab(NULL) +
#      ylab(NULL) +
#      scale_fill_manual(values = c("#050505", "#E69F00", "#56B4E9"),
#  		name = "Sample Size") +
#      facet_wrap(factor(truedf)~.,  scales="free") +
#      NULL
#  denseplot
#  
#  
#  knitr::kable(dfdataframe %>% group_by(truedf, samplesize) %>%
#  	summarize(min = min(estdf), max = max(estdf),
#  	median = median(estdf),
#  	mean=mean(estdf),
#  	sd = sd(estdf)))
#  ### Here ends the long simulation
#  
#  ##### Here is what is really run
#  
#  set.seed(20190621)
#  df10 <- rep(0,50)
#  df1050 <-  rep(0,50)
#  df10100 <-  rep(0,50)
#  
#  for(i in 1:50){
#     df10[i] = suppressWarnings(MLmatrixt(rmatrixt(mean = matrix(0,5,3),df = 10, n = 25), fixed = FALSE, df = 5, max.iter = 20)$nu)
#     df1050[i] = suppressWarnings(MLmatrixt(rmatrixt(mean = matrix(0,5,3),df = 10, n = 50), fixed = FALSE, df = 5, max.iter = 20)$nu)
#     df10100[i] = suppressWarnings(MLmatrixt(rmatrixt(mean = matrix(0,5,3),df = 10, n = 100), fixed = FALSE, df = 5, max.iter = 20)$nu)
#  }
#  
#  
#  dfdataframe = data.frame(label = c('10 df'),
#                           estdf = c(df10, df1050, df10100),
#                           samplesize = factor(rep(c(25,50,100), each = 50)))
#  library(ggplot2)
#  library(dplyr)
#  library(magrittr)
#  denseplot <- ggplot(data = subset(dfdataframe, estdf < 200),aes(x=estdf, fill=samplesize)) +
#      geom_density(alpha = .5) +
#      geom_vline(mapping = aes(xintercept = 10),
#                 size = .5) +
#      theme_bw() +
#      theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), strip.text = element_text(size = 8),
#            legend.justification=c(1,0), legend.position=c(.95,.4),
#            legend.background = element_blank(),
#            legend.text =element_text(size = 8), legend.title = element_text(size = 8)) +
#      ggtitle("Density plot of estimated degrees of freedom compared to actual") +
#      xlab(NULL) +
#      ylab(NULL) +
#      scale_fill_manual(values = c("#050505", "#E69F00", "#56B4E9"), name = "Sample Size") +
#   #   facet_wrap(factor(truedf)~.,  scales="free") +
#      NULL
#  denseplot
#  
#  
#  knitr::kable(dfdataframe %>% group_by(samplesize) %>%
#  	summarize(min = min(estdf), max = max(estdf),
#  	median = median(estdf),
#  	mean=mean(estdf),
#  	sd = sd(estdf)))
#  
#  #### Here ends what is really run
#  
#  A <- rmatrixt(30, mean = matrix(0, nrow=2, ncol=3), df = 10)
#  B <- rmatrixt(30, mean = matrix(c(1,0), nrow=2, ncol=3), df = 10)
#  C <- rmatrixt(30, mean = matrix(c(0,1), nrow=2, ncol=3), df = 10)
#  ABC <- array(c(A,B,C), dim = c(2,3,90))
#  groups <- factor(c(rep("A",30),rep("B",30),rep("C",30)))
#  prior = c(30,30,30)/90
#  matlda <- matrixlda(x = ABC,grouping = groups, prior = prior,
#                      method = 't', nu = 10, fixed = TRUE)
#  predict(matlda, newdata = ABC[,,c(1,31,61)])
#  
#  sessionInfo()

