\name{groupedSurv}
\alias{groupedSurv}
\title{
    Compute the Efficient Score Statistic for the Grouped Survival Model
}
\description{
  A method to compute the efficient score statistic for a grouped survival model. The function supports zero, single, or multiple covariates as dictated by the input data. The function can take either a \code{gwaa.data} object or a standard \code{data.frame} or \code{matrix} as input. If more than one variable of interest is provided, the efficient score will be computed for each one separately. For each variable tested, the function returns the unadjusted marginal asymptotic p-value, as well as the family-wise error rate (FWER, Bonferroni correction) adjusted p-value and local false-discovery rates (FDR, Storey's Q-values).
}

\usage{
  groupedSurv(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL, 
    gtime, delta, beta=0, nCores=1, reScore=FALSE)
}

\arguments{
  \item{x}{Vector, \code{data.frame}, or \code{matrix} of numeric variables of interest for each sample. Alternatively, if a \code{gwaa.data} object is given for \code{GenABEL.data}, this argument may be a vector of strings corresponding to column names in \code{GenABEL.data@gtdata} to use as variables of interest, set \code{NULL} will use all the variables in \code{GenABEL.data@gtdata} as variables of interest.}  
  \item{Z}{Optional \code{data.frame} or \code{matrix} of numeric covariate values for each sample. Alternatively, if a \code{gwaa.data} object is given for \code{GenABEL.data}, this argument may be a vector of strings corresponding to column names in \code{GenABEL.data@phdata} to use as covariates.}
  \item{GenABEL.data}{Optional \code{gwaa.data} object.}
  \item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{theta}{Optional vector of estimated nuisance parameters for the covariates.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
  \item{beta}{Scalar for the parameter of interest. Default is 0.}
  \item{nCores}{Integer representing the number of cores to be used for multi-threaded computation. Default is 1.}
  \item{reScore}{Boolean indicating whether to return the full efficient scores matrix. Default is \code{FALSE}}
}


\value{
  If \code{reScore=FALSE}, the function returns a \code{data.frame} with one row for each variables of interest and four columns: \code{stat}, the efficient score statistic, \code{pvalue}, unadjusted p-value, \code{FWER}, the family-wise error rate adjusted p-value (Bonferroni correction), and \code{FDR}, the local false-discovery rate (Storey's Q-value). If \code{reScore=TRUE}, a \code{list} is returned. The first element is the \code{data.frame} of efficient score statistics and p-values, and the second is a matrix of contributions of each sample to the efficient score statistics of each variable of interest. Rows correspond to samples, and columns correspond to the variables of interest.
}

\references{
  Bonferroni, C.E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in Onore del Professore Salvatore Ortu Carbon, 13-60.\cr
  Storey, J.D. (2003). The positive false discovery rate: a Bayesian interpretation and the q-value. Annals of Statistics, 31:6, 2013-2035.\cr
  Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: a unified approach. Journal of the Royal Statistical Society Series B-Statistical Methodology, 66, 187-205\cr
}


\examples{
# Generate dummy data
set.seed(111)
n <- 1000

# effect size
beta <- 0.4

# covariate parameters
theta <- c(0.2, 0.2)

# variable of interest associated with outcome
MAF <- 0.05
x <- matrix(rbinom(n, 2, MAF), ncol = 1)

# covariate data (centered at 0)
z1 <- rnorm(n)
z2 <- rbinom(n, 1, 0.5) - 0.5
Z  <- matrix(cbind(z1, z2), ncol = 2)

# continuous survival time
lam0 <- 1
cmax <- 3
lami <- lam0 * exp(x*beta + Z[,1]*theta[1]+Z[,2]*theta[2])
stime <- rexp(n, lami)
ctime <- runif(n, 0, cmax)
delta <- stime < ctime
otime <- pmin(stime, ctime)

# additional variables of interest
xMore <- matrix(rbinom(n*100, 2, MAF), ncol = 100)
xMore <- cbind(x, xMore)

# number of observation times
ntps <- 5

# number of intervals
r <- ntps + 1

# last observation time
maxbreakq <- 0.85
maxbreak <- qexp(maxbreakq, lam0)

# grouped failure times
breaks <- (1:ntps) * (maxbreak/ntps)
gtime <- findInterval(otime, breaks) + 1
delta[gtime == r] <- FALSE
dctime <- findInterval(ctime, breaks) + 1
delta[gtime == dctime] <- FALSE
delta <- as.numeric(delta)
gtime[which(gtime == r)] <- Inf

# estimate nuisance parameters
#thetaest <- thetaEst(Z, gtime, delta)

# compute efficient score statistics
# eff <- groupedSurv(x=xMore, Z=Z, alpha=thetaest$alpha, theta=thetaest$theta,
#		gtime=gtime, delta=delta, beta=0, nCores=1)
# head(eff)
}

