## ----setup1, include=FALSE, echo=FALSE----------------------------------------
require(knitr)

## ----setup2, include=FALSE,echo=FALSE-----------------------------------------
library(groupedSurv)
old_ops <- options(width=80)  # make the printing fit on the page
set.seed(1121)     # make the results repeatable
stdt<-date()

## ----without, eval=FALSE------------------------------------------------------
#  thetaEst(Z=NULL, gtime, delta, method="BFGS")
#  betaEst(x, Z=NULL, alpha, theta=NULL, gtime, delta)
#  groupedSurv(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL,
#           gtime, delta, beta=0, nCores=1, reScore=FALSE)
#  geneStat(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL,
#           gtime, delta, beta=0, nCores=1,
#           FUN=function(Uij, weight){sum((colSums(Uij)*weight)^2)}, geneSet)

## ----with, eval=FALSE---------------------------------------------------------
#  alphaEstFam(gtime, delta)
#  betaEstFam(x, fam_group, fam_role, alpha, var, gtime, delta, lower, upper)
#  varEstFam(x, fam_group, fam_role, alpha, gtime, delta, lower, upper, beta = 0)
#  groupedSurvFam(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0)
#  varLLFam(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0)

## ----inputParawithoutEffScore-------------------------------------------------
set.seed(111)
n <- 1000
# effect size
beta <- 0.3
# covariate parameters
theta <- c(0.2, 0.2) 
# variable of interest associated with outcome
MAF  <- 0.05
x <- matrix(rbinom(n, 2, MAF), ncol = 1)
# additional variables of interest
xMore <- matrix(rbinom(n*100, 2, MAF), ncol = 100)
xMore <- cbind(x, xMore)
# covariate data (centered at 0) 
z1 <- rnorm(n)
z2 <- rbinom(n, 1, 0.5) - 0.5
Z <- matrix(cbind(z1, z2), ncol = 2)
# continuous survival time
lam0 <- 1
cmax <- 3
lami <- lam0 * exp(x * beta + Z %*% theta)
stime <- rexp(n, lami)
ctime <- runif(n, 0, cmax)
delta <- stime < ctime
otime <- pmin(stime, ctime)

## ----inputParawithoutEffScore2,eval=TRUE--------------------------------------
# grouped observation time points
ntps <- 5
r <- ntps + 1
# last observation time
maxbreakq <- 0.85
maxbreak <- qexp(maxbreakq, lam0)
# grouped survival times
breaks <- (1:ntps) * (maxbreak/ntps)
gtime <- findInterval(otime, breaks) + 1 
delta[gtime == r] <- FALSE 
dctime <- findInterval(ctime, breaks) + 1
delta[gtime == dctime] <- FALSE
delta <- as.numeric(delta)
gtime[which(gtime == r)] <- Inf
table(gtime, delta)

## ----loadpkg,eval=TRUE--------------------------------------------------------
library(groupedSurv)

## ----expa1,eval=FALSE---------------------------------------------------------
#  thetaest <- thetaEst(Z, gtime, delta)
#  thetaest

## ----expa2, eval=FALSE--------------------------------------------------------
#  eff <- groupedSurv(x=xMore, Z=Z,  alpha=thetaest$alpha, theta=thetaest$theta,
#                     gtime=gtime, delta=delta, beta=0, nCores=1)
#  head(eff)

## ----expa4, eval=FALSE--------------------------------------------------------
#  betaest <- betaEst(x=x, Z=Z, alpha=thetaest$alpha, theta=thetaest$theta,
#                     gtime=gtime, delta=delta)
#  betaest

## ----expa3, eval=FALSE, message=FALSE-----------------------------------------
#  library(GenABEL)
#  data(srdta)
#  GenABELdat <- srdta[1:n]
#  snpsToTest <- GenABELdat@gtdata@snpnames[1:200]
#  eff <- groupedSurv(x=snpsToTest, Z=Z, GenABEL.data=GenABELdat,
#                     alpha=thetaest$alpha, theta=thetaest$theta,
#                     gtime=gtime, delta=delta, beta=0, nCores=1)

## ----bedtoSNPmatrix, eval=FALSE-----------------------------------------------
#  library(BEDMatrix)
#  path <- system.file("extdata", "example.bed", package = "BEDMatrix")
#  m <- BEDMatrix(path)
#  # Extract genotypes for the specified variants
#  xPLINK <- m[, c("snp0_A", "snp1_C", "snp2_G")]

## ----VCFtoSNPmatrix, eval=FALSE-----------------------------------------------
#  system("wget ftp://share.sph.umich.edu/minimac3/DosageConvertor/DosageConvertor.v1.0.4.tar.gz")
#  system("tar -xzvf DosageConvertor.v1.0.4.tar.gz")
#  library(VariantAnnotation)
#  exampleVcfFile <- "./DosageConvertor/test/TestDataImputedVCF.dose.vcf.gz"
#  myvcf <- readVcf(exampleVcfFile, "hg19")
#  dosedat <- assay(myvcf,"DS")
#  xVCF <- t(dosedat)

## ----geneStat3, eval=TRUE-----------------------------------------------------
geneInfo <- data.frame(gene=c("BRCA1","BRCA2"), chr=c(17,13),
                       start=c(41196312, 32889611), end=c(41277500, 32973805), 
                       stringsAsFactors=FALSE)
snpInfo <- data.frame(chr=c(17,17,13,13), pos=c(41211653,41213996,32890026,32890572),
                      rsid=c("rs8176273","rs8176265","rs9562605","rs1799943"),
                      stringsAsFactors=FALSE)

## ----snplist, eval=FALSE,  results='hide'-------------------------------------
#  library(snplist)
#  setGeneTable(geneInfo)
#  setSNPTable(snpInfo)
#  geneset <- makeGeneSet()

## ----geneset, eval=FALSE------------------------------------------------------
#  geneset

## ----geneStat4, eval=FALSE----------------------------------------------------
#  G <- matrix(rbinom(n*nrow(snpInfo), 2, 0.5), ncol=nrow(snpInfo))
#  colnames(G) <- snpInfo$rsid

## ----geneStat5, eval=FALSE----------------------------------------------------
#  for(i in seq_len(length(geneset))){
#    weight <- rep(1, length(geneset[[i]]))
#    geneset[[i]] <- list(geneset[[i]], weight)
#  }

## ----geneStat, eval=FALSE-----------------------------------------------------
#  res <- geneStat(x=G, Z=Z, alpha=thetaest$alpha, theta=thetaest$theta,
#                  gtime=gtime, delta=delta, geneSet=geneset)
#  res$stat

## ----inputPara3---------------------------------------------------------------
#rm(list=ls())
set.seed(111)
m <- 10

# family ID
fgrp <- as.character(rep(1:m, each=3))

# role within family
f_ind <- rep(c('o','f','m'),m)

# grouped survival data
gtimes <- sample(1:4, m*3, replace=TRUE)
deltas <- sample(0:1, m*3, replace=TRUE)

# variable of interest
g <- rbinom(m*3, 2, 0.1)

# parameter search bounds
upper  <- 2
lower  <- 0

## ----expa, eval=FALSE---------------------------------------------------------
#  alphaest <- alphaEstFam(gtimes, deltas)
#  alphaest

## ----expv, eval=FALSE---------------------------------------------------------
#  varest <- varEstFam(x=g, fam_group=fgrp, fam_role=f_ind, alpha=alphaest,
#                      gtime=gtimes, delta=deltas, lower, upper, beta=0)
#  varest

## ----effscore3, eval=FALSE----------------------------------------------------
#  effFam <- groupedSurvFam(x=g, fam_group=fgrp, fam_role=f_ind, alpha=alphaest,
#                        var=varest, gtime=gtimes, delta=deltas, beta=0)
#  PvalueFam(effFam)

## ----expb2, eval=FALSE--------------------------------------------------------
#  betaEstFam(x=g, fam_group=fgrp, fam_role=f_ind, alpha=alphaest,
#             var=varest, gtime=gtimes, delta=deltas, lower, upper)

## ----sessinfo, echo=FALSE, include=TRUE, results='asis'-----------------------
toLatex(sessionInfo(), locale=FALSE)
## reset options
options(old_ops)

