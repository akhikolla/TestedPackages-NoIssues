\name{geneStat}
\alias{geneStat}
\title{
  Compute Gene-Level Statistics Based on Efficient Score Statistics for the Grouped Survival Model
}
\description{
  A method to compute gene-level statistics using the contributions of individual samples to the efficient score statistics of SNPs under a grouped survival model. The function supports zero, single, or multiple covariates as dictated by the input data. The function can take either a \code{gwaa.data} object or a standard \code{data.frame} or \code{matrix} as input. The contribution of each sample to the efficient score will be computed separately for each SNP being tested, and then gene-level statistics will be computed using either the default Sequence Kernel Association Test (SKAT) or a user-provided function.
}

\usage{
	geneStat(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL, gtime, delta,
           beta=0, nCores=1, 
           FUN=function(Uij,weight){sum((colSums(Uij)*weight)^2)}, geneSet)
}

\arguments{
  \item{x}{Vector, \code{data.frame}, or \code{matrix}  of numeric variables of interest for each sample. Alternatively, if a \code{gwaa.data} object is given for \code{GenABEL.data}, this argument may be a vector of strings corresponding to column names in \code{GenABEL.data@gtdata} to use as variables of interest.}  
  \item{Z}{Optional \code{data.frame} or \code{matrix} of numeric covariate values for each sample. Alternatively, if a \code{gwaa.data} object is given for \code{GenABEL.data}, this argument may be a vector of strings corresponding to column names in \code{GenABEL.data@phdata} to use as covariates.}
  \item{GenABEL.data}{Optional \code{gwaa.data} object.}
  \item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{theta}{Optional vector of estimated nuisance parameters for the covariates.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
  \item{beta}{Scalar for the parameter of interest. Default is 0.}
  \item{nCores}{Integer representing the number of cores to be used for multi-threaded computation. Default is 1.}
  \item{FUN}{A function to compute the gene-level statistics base on the SNP-level efficient scores. The function arguments should include a matrix of each samples' contribution to the efficient scores, \code{Uij}, and a vector of weights, \code{weight}. Default is the SKAT statistic.}
  \item{geneSet}{A list of lists, where each sub-list contains two vectors: a vector of column names of the variables of interest to be included in that set, and a vector of numerical weights corresponding to each specified variable of interest.}
}

\value{
  \code{list} with two elements: \code{stat}, a list of return objects from the default or user-defined \code{FUN}, and \code{U}, a list of matrices of the contribution of each subject to the efficient score statistics for each SNP in each gene set. Rows correspond to samples, and columns correspond to the variables of interest.
}

\references{
	Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M. and Lin, X. (2011). Rare-variant association testing for sequencing data with the sequence kernel association test. American Journal of Human Genetics, 89:1, 82-93.\cr
  Ionita-Laza, I., Lee, S., Makarov, V., Buxbaum, J.D., and Lin, X. (2013). Sequence kernel association tests for the combined effect of rare and common variants. American Journal of Human Genetics, 92:6, 841-53.\cr
}

\examples{
# Generate dummy data
rm(list=ls())
set.seed(111)
n <- 1000

# covariate parameters
theta <- c(0.2, 0.2)

# covariate data (centered at 0)
z1 <- rnorm(n)
z2 <- rbinom(n, 1, 0.5) - 0.5
Z <- matrix(cbind(z1, z2), ncol = 2)

# continuous survival time
lam0 <- 1
cmax <- 3
lami <- lam0 * exp(Z[,1]*theta[1]+Z[,2]*theta[2])
stime <- rexp(n, lami)
ctime <- runif(n, 0, cmax)
delta <- stime < ctime
otime <- pmin(stime, ctime)

# number of observation times
ntps <- 5

# number of intervals
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

# estimate nuisance parameters
#thetaest <- thetaEst(Z, gtime, delta)

# SNP and gene information
#geneInfo <- data.frame(gene=c("BRCA1","BRCA2"), chr=c(17,13),
#                       start=c(41196312, 32889611), end=c(41277500, 32973805),
#                       stringsAsFactors=FALSE)
#snpInfo <- data.frame(chr=c(17,17,13,13),
#                      pos=c(41211653,41213996,32890026,32890572),
#                      rsid=c("rs8176273","rs8176265","rs9562605","rs1799943"),
#                      stringsAsFactors=FALSE)

# use snplist package to create gene sets
#library(snplist)
#setGeneTable(geneInfo)
#setSNPTable(snpInfo)
#geneset <- makeGeneSet()

# simulate sample values for SNPs
#G <- matrix(rbinom(n*nrow(snpInfo), 2, 0.5), ncol=nrow(snpInfo))
#colnames(G) <- snpInfo$rsid

# dummy weights for the SNPs in the gene sets
#for(i in seq_len(length(geneset))){
#    weight <- rep(1, length(geneset[[i]]))
#    geneset[[i]] <- list(geneset[[i]], weight)

#}

# compute gene-level statistics based on the default SKAT function
#res <- geneStat(x=G, Z=Z, alpha=thetaest$alpha, theta=thetaest$theta,
#                gtime=gtime, delta=delta, geneSet=geneset, beta=0, nCores=1)
#res$stat
}










