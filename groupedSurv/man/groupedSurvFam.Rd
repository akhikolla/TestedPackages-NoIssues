\name{groupedSurvFam}
\alias{groupedSurvFam}
\title{
	Compute the Efficient Score Statistic for the Frailty Model
}
\description{
	A method to compute the efficient score statistic for a frailty model,
	accounting for family structure of related individuals (e.g., trios). The
	input data is assumed to be organized such that records for each family occur
	consecutively, and that records for offspring precede those for parents. The
	variance matrix for the random effects is assumed to be of the form
	\code{var*K}, where \code{K} is a matrix of kinship coefficients between
	family members. The following groupings are permitted: (Individual),
	(Offspring, Offspring), (Offspring, Parent), (Offspring, Parent, Parent), and
	(Offspring, Offspring, Parent, Parent). Other family structures have not been
	implemented. 
}

\usage{
  groupedSurvFam(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0, nCores=1)  	
}

\arguments{
  \item{x}{Vector or \code{matrix} of numeric variables of interest for each sample.}
  \item{fam_group}{Vector of family IDs for each sample.}
  \item{fam_role}{Vector of indicators for the role within a family of each sample, i.e., \{"Offspring","Mother","Father"\}, or \{"o","m","f"\}.}
  \item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{var}{Scalar for frailty variance.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
  \item{beta}{Scalar for the fixed effect size. Default is 0.}
  \item{nCores}{Integer representing the number of cores to be used for multi-threaded computation. Default is 1.}
}

\value{
	Matrix of contributions of each family to the efficient score statistics of each variable of interest. Rows correspond to families, and columns correspond to the variables of interest.
}

\examples{
# Generate dummy data
set.seed(111)
n <- 24

# effect size
beta <- 0.3

# covariate parameters
theta <- c(0.2, 0.2)

# variable of interest
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

# family-specific information
m <- n/3
fam_role <- rep(c("o","f","m"),m)
fam_group <- as.character(rep(1:m,rep(3,m)))

# nuisance parameters
#alpha <- thetaEst(Z, gtime,delta)$alpha
#var <- 0.2

# compute efficient score statistics
#res   <- groupedSurvFam(x, fam_group, fam_role, alpha, var, gtime, delta) 
#res
}

