## ----setup, include=FALSE-----------------------------------------------------
# library(rgl)
# #library(rglwidget)
# setupKnitr()
# knitr::opts_chunk$set(echo = TRUE,
#                       fig.align = "center",
#                       warning = FALSE,
#                       webgl = TRUE,
#                       fig.width = 8, 
#                       fig.height = 8,
#                       fig.keep = "all",
#                       fig.ext = "jpeg"
#                       )

## ----fig.width=5, fig.height=5------------------------------------------------
library(DataVisualizations)
data("Lsun3D")
Pixelmatrix(Lsun3D$Data)

## -----------------------------------------------------------------------------
library(DataVisualizations)
data(MTY)
InspectVariable(MTY,'MTY')


## ----fig.width=4, fig.height=4, message=FALSE---------------------------------
library(DataVisualizations)
library(ggplot2)
data(ITS)
data(MTY)
library(vioplot)
Data=cbind(ITS,MTY)
#MDplot(Data)+ylim(0,6000)+ggtitle('Two Features With Adjusted Range')

#MDplot(Data,Scaling = "Robust")+ggtitle('"Shape-Invariant" Normalization')

#Data is now capped
#Data[Data[,2]>6000,2]=6000
MDplot(Data)+ylim(0,6000)+ggtitle('Two Features with MTY Capped')
boxplot(Data,main='Two Features with MTY Capped')
vioplot(Data[,1],Data[,2])
title('Two Features with MTY Capped')

## -----------------------------------------------------------------------------
library(DataVisualizations)
data(ITS)
data(MTY)
Ind2=which(ITS<900&MTY<8000)
PDEscatter(ITS[Ind2],MTY[Ind2],xlab = 'ITS in EUR',ylab ='MTY in EUR' ,main='Scatter density plot using PDE' )

## ----fig.width=4, fig.height=4------------------------------------------------
data("Lsun3D")
n=nrow(Lsun3D$Data)
Data=cbind(Lsun3D$Data,runif(n),rnorm(n),rt(n,2),rlnorm(n),rchisq(100,2))
Header=c('x','y','z','uniform','gauss','t','log-normal','chi')
cc=cor(Data,method='spearman')
diag(cc)=0
Pixelmatrix(cc,YNames = Header,XNames = Header,main = 'Spearman Coeffs')

## -----------------------------------------------------------------------------
library(DataVisualizations)
data("Lsun3D")
InspectDistances(Lsun3D$Data,method="euclidean")

## -----------------------------------------------------------------------------
library(DataVisualizations)
Cls=c(1L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 
2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 
1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 
2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 3L, 2L, 2L, 2L, 1L, 
2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 
2L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 
2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 
2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 
2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L
)
Codes=c("AFG", "AGO", "ALB", "ARG", "ATG", "AUS", "AUT", "BDI", "BEL", 
"BEN", "BFA", "BGD", "BGR", "BHR", "BHS", "BLZ", "BMU", "BOL", 
"BRA", "BRB", "BRN", "BTN", "BWA", "CAF", "CAN", "CH2", "CHE", 
"CHL", "CHN", "CIV", "CMR", "COG", "COL", "COM", "CPV", "CRI", 
"CUB", "CYP", "DJI", "DMA", "DNK", "DOM", "DZA", "ECU", "EGY", 
"ESP", "ETH", "FIN", "FJI", "FRA", "FSM", "GAB", "GBR", "GER", 
"GHA", "GIN", "GMB", "GNB", "GNQ", "GRC", "GRD", "GTM", "GUY", 
"HKG", "HND", "HTI", "HUN", "IDN", "IND", "IRL", "IRN", "IRQ", 
"ISL", "ISR", "ITA", "JAM", "JOR", "JPN", "KEN", "KHM", "KIR", 
"KNA", "KOR", "LAO", "LBN", "LBR", "LCA", "LKA", "LSO", "LUX", 
"MAC", "MAR", "MDG", "MDV", "MEX", "MHL", "MLI", "MLT", "MNG", 
"MOZ", "MRT", "MUS", "MWI", "MYS", "NAM", "NER", "NGA", "NIC", 
"NLD", "NOR", "NPL", "NZL", "OMN", "PAK", "PAN", "PER", "PHL", 
"PLW", "PNG", "POL", "PRI", "PRT", "PRY", "ROM", "RWA", "SDN", 
"SEN", "SGP", "SLB", "SLE", "SLV", "SOM", "STP", "SUR", "SWE", 
"SWZ", "SYC", "SYR", "TCD", "TGO", "THA", "TON", "TTO", "TUN", 
"TUR", "TWN", "TZA", "UGA", "URY", "USA", "VCT", "VEN", "VNM", 
"VUT", "WSM", "ZAF", "ZAR", "ZMB", "ZWE")
Worldmap(Codes,Cls)

## ----fig.width=5, fig.height=5------------------------------------------------
library(DataVisualizations)
data(categoricalVariable)
Fanplot(categoricalVariable)

Piechart(categoricalVariable)

## ----warning=FALSE, comment=FALSE---------------------------------------------
library(DataVisualizations)
data("Lsun3D")

Heatmap(Lsun3D$Data,Lsun3D$Cls,method = 'euclidean')

Silhouetteplot(Lsun3D$Data,Lsun3D$Cls,PlotIt = T)

## ----fig.width=4, fig.height=4,warning=FALSE----------------------------------
library(DataVisualizations)
data("Lsun3D")
Accuracy=matrix(NaN,100,2)
Algorithms=c("MacQueen","Lloyd")
colnames(Accuracy)=Algorithms
for(i in 1:100){
  Cls=kmeans(Lsun3D$Data,4,algorithm=Algorithms[1])$cluster
  Cls2=kmeans(Lsun3D$Data,4,algorithm=Algorithms[2])$cluster
  #this is an artifical example, because the problem of arbitrary class labels is not accounted for
  #please choose an appropiate internal index or an external index
  Accuracy[i,1]=sum(Cls==Lsun3D$Cls)/length(Lsun3D$Cls)
  Accuracy[i,2]=sum(Cls2==Lsun3D$Cls)/length(Lsun3D$Cls)
}

MDplot(Accuracy) + xlab('Output of Evaluation of two Algorithms') +
  ylab('Range of Values of the Evaluation of an Algorithm') +
  ggtitle("Simple Benchmarking")

