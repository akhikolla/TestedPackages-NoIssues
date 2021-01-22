\name{covglasso}
\alias{covglasso}
\title{Sparse covariance matrix estimation}

\description{
Direct estimation of a sparse covariance matrix using the covariance graphical lasso.
}

\usage{
covglasso(data = NULL,
          S = NULL, n = NULL,
          lambda = NULL,
          rho = NULL,
          duplicated = TRUE,
          L = 10,
          crit = c("bic", "ebic"),
          gamma = 0.8,
          penalize.diag = FALSE,
          start = NULL,
          ctrl = control(),
          path = FALSE)
}

\arguments{
\item{data}{A numerical dataframe or matrix, where rows correspond to observations and columns to variables. If \code{data = NULL}, the sample covariance \code{S} must be provided in input.}

\item{S}{The sample covariance matrix of the data. If \code{S = NULL}, the maximum likelihood estimate of the covariance matrix is used in the estimation of the sparse covariance matrix.}

\item{n}{The number of observations. If \code{data = NULL} and \code{S} is provided in input, \code{n} must be provided in input as well.}

\item{lambda}{A vector or array of non-negative lasso regularization parameters. Penalization is applied elementwise to all entries of the covariance matrix. If an array, each entry must be a matrix with same dimensions of the sample covariance matrix. Values should be increasing from the smallest to the largest. If \code{lambda = NULL}, an alternative penalization based on thresholding of the empirical correlation matrix is used; see "Details".}

\item{rho}{A vector of correlation values used to define the penalization in terms of the thresholded sample correlation matrix. See "Details". Note that this penalization is used by default.}

\item{duplicated}{Remove duplicated penalty matrices when the default penalty term based on the thresholded correlation matrix is used. Suggest to leave this argument to \code{TRUE} all the time as several redundant matrices giving the same penalty term are discarded.}

\item{L}{The number of \code{rho} values. Only used when \code{lambda} and \code{rho} are NULL. Default is \code{L = 10}.}

\item{crit}{The model selection criterion employed to select the optimal covariance graph model. Can be \code{"bic"} or \code{"ebic"}; see "Details".}

\item{gamma}{A penalty parameter used when \code{crit = "ebic"} and EBIC is used to select the optimal graph covariance model. The value of \code{gamma} must be in the range \code{[0,1]}. Default is \code{gamma = 0.8}, which encourages sparser models.}

\item{penalize.diag}{A logical argument indicating if the diagonal of the covariance matrix should be penalized. Default to \code{FALSE}.}

\item{start}{A starting matrix for the estimation algorithm. If \code{NULL}, the starting value is the diagonal sample covariance matrix.}

\item{ctrl}{A list of control parameters for the coordinate descent algorithm employed for estimation. See also \code{\link{control}}.}

\item{path}{A logical argument controlling whether all the estimated covariance matrices along the path defined by \code{lambda} or \code{rho} should be included in the output.}
}

\value{
A list containing the following elements.
\item{sigma}{The estimated covariance matrix.}
\item{omega}{The estimated concentration (inverse covariance) matrix.}
\item{graph}{The adjacency matrix given in input corresponding to the marginal or conditional independence graph.}
\item{loglik}{Value of the maximized log-likelihood.}
\item{npar}{Number of estimated non-zero parameters.}
\item{penalty}{Value of the penalty term.}
\item{bic}{Optimal BIC or EBIC value.}
\item{BIC}{All BIC or EBIC values along the path defined by \code{lambda} or \code{rho}.}
\item{path}{A list containing all the estimated sparse covariance models. Provided in output only when \code{path = TRUE}.}
\item{rho}{The values of \code{rho} thresholds used to define the penalization based on the thresholded sample correlation matrix.}
\item{lambda}{The values of \code{lambda} penalty parameters for the penalization.}
}

\details{
The function estimates a sparse covariance matrix using a fast coordinate descent algorithm to solve the covariance graphical lasso. The estimated sparse covariance matrix is obtained by optimizing the following penalized log-likelihood:
\deqn{-\frac{n}{2}\{ \mathrm{logdet}(\Sigma) + \mathrm{trace}(S\Sigma^{-1})  \} - ||\Lambda * \Sigma||_1}{-n/2*( logdet(Sigma) + trace(S \%*\% Sigma^-1) ) - ||lambda*Sigma||_1}
subject to \eqn{\Sigma}{Sigma} being positive definite. In the penalty term, the \eqn{L_1}{L1} norm and the matrix multiplication between \eqn{\Lambda}{lambda} and \eqn{\Sigma}{Sigma} is elementwise.

By default (when \code{lambda = NULL}), the penalization matrix \eqn{\Lambda}{lambda} is defined in terms of a sequential thresholding of the sample correlation matrix. Given \eqn{\rho_l}{rho[l]} a threshold value and \eqn{R}{R} the sample correlation matrix, the penalty term matrix \eqn{\Lambda}{lambda} is defined by the values \eqn{(1/s_{ij})I(r_{ij} < \rho_l)}{ifelse(r[i,j] < rho[l], 1/s[i,j], 0)}, that is:
\deqn{\Lambda = \frac{1}{S}I(R < \rho_l)}{lambda = ifelse(R < rho[l], 1/S, 0)}
where the inequality is taken elemntwise. Such choice of penalty matrix provides a framework related to the method of Chaudhuri et al. (2007).
If the vector \code{rho} is not given in input, the sequence of threshold values is defined as the \code{L} quantiles of the absolute values of the sample correlations in \eqn{R}{R}. If \code{lambda} is provided in input, the penalization corresponds to the standard covariance graphical lasso of Bien, Tibshirani (2011).

The sparse covariance matrix corresponds to a Gaussian covariance graphical model of marginal independence, where in the sparse covariance matrix a zero entry corresponds to two variables being marginally independent. Different penalizations \code{lambda} imply different models, and selection of the optimal graphical model is performed using \code{"bic"} (default) or \code{"ebic"}. In the latter case, the argument \code{gamma} controls the additional penalty term in the model selection criterion; see Foygel, Drton, (2010).
}


\examples{
# a simple example with a 3-block diagonal matrix
library(MASS)
v <- 3
n <- 300
sig <- matrix(0.8, v,v)
diag(sig) <- 1
set.seed(190188)
tmp <- replicate( 3, mvrnorm(n, rep(0,v), sig) )
x <- matrix(c(tmp), n, v*3)

fit1 <- covglasso(x)
plot(fit1$rho, fit1$BIC)
image(fit1$sigma != 0)

# refine search
fit2 <- covglasso(x, rho = seq(0.1, 0.4, length = 50))
image(fit2$sigma != 0)

fit1$bic
fit2$bic


# Cars93 data in MASS package
data("Cars93", package = "MASS")
dat <- na.omit( Cars93[,c(4:8,12:15,17,19:25)] )

fit1 <- covglasso(dat, L = 50)

# more sparse
fit2 <- covglasso(dat, L = 50,
                    crit = "ebic", gamma = 1)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1,2))
plot(fit1$rho, fit1$BIC, main = "BIC")
plot(fit2$rho, fit2$BIC, main = "EBIC")
image(fit1$sigma != 0, col = c("white", "black"), main = "BIC")
image(fit2$sigma != 0, col = c("white", "black"), main = "EBIC")
par(oldpar)  # reset par

\donttest{
# a more complex situation
library(mixggm)

V <- 20
M <- choose(V, 2)
n <- 200
rand <- rbinom(M, 1, 0.4)
graph <- matrix(0, V,V)
graph[upper.tri(graph)] <- rand
graph <- graph + t(graph)
sigma <- matrix(0.8, V,V)
diag(sigma) <- 1
sigma <- fitGGM(data = NULL, S = sigma, n, graph = graph)$sigma
x <- mvrnorm(n, rep(0, V), sigma)

fit <- covglasso(x, L = 50)
plot(fit$rho, fit$BIC)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1,2))
cols <- c("white", grey.colors(10, start = 0, end = 0.7, rev = TRUE))
image(fit$sigma, main = "Fitted", col = cols)
image(sigma, main = "Actual", col = cols)
par(oldpar)   # reset par
}
}


\references{
Bien, J., Tibshirani, R.J. (2011). Sparse estimation of a covariance matrix. \emph{Biometrika}, 98(4), 807–820.

Chaudhuri, S., Drton M., Richardson, T. S. (2007). Estimation of a covariance matrix with zeros. \emph{Biometrika}, 94(1), 199-216.

Foygel, R., Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. In \emph{Advances in neural information processing systems}, pages 604–612.

Wang, H. (2014). Coordinate descent algorithm for covariance graphical lasso. \emph{Statistics and Computing}, 24:521.
}


\seealso{
\code{\link{control}}
}
