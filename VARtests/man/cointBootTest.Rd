\name{cointBootTest}

\alias{cointBootTest}
\alias{print.cointBootTest}

\title{
Bootstrap Determination of Cointegration Rank in VAR Models
}

\encoding{UTF-8}

\description{
This function uses the bootstrap and wild bootstrap to test the cointegration rank of a VAR model. 
The test is an implementation of Cavaliere, Rahbek & Taylor (2012, 2014), and is used in Ahlgren & Catani (2018).
}
\usage{
cointBootTest(y, r = "sequence", p, model = 1, signif = 0.05, dummies = NULL, B = 999,
boot_type = c("B", "WB"), WB_dist = c("rademacher", "normal", "mammen"))
\method{print}{cointBootTest}(x, ...)
}



\arguments{
  \item{y}{
a T x K matrix containing the time series.
}
  \item{r}{
either \code{"sequence"} or a vector of integers (0 <= r <= K - 1, where K is the number of columns in \code{y}).
If a vector of integers, \code{r} is the cointegration rank being tested.
If \code{r = "sequence"}, the bootstrap sequential algorithm will be used (see 'details').
}
  \item{p}{
the lag order of the model.
}
  \item{model}{
either 1 (no deterministic terms), 2 (restricted constant), or 3 (restricted linear trend). See 'details' below.
}
  \item{signif}{
if \code{r = "sequence"} is used, \code{signif} sets the significance level of the tests in the sequential algorithm.
}
  \item{dummies}{
(optional) dummy variables. Must have the same number of rows as \code{y}. The models will then be estimated with the dummy variables, 
but the dummy variables are not used in the bootstrap DGP. In the \code{boot_type = "B"} version, the residuals used to draw the bootstrap errors
do not include rows corresponding to observations where any of the dummies are equal to 1.
}
  \item{B}{
the number of bootstrap replications.
}
  \item{boot_type}{
either "B", "WB", or both. "B" uses the iid bootstrap algorithm, while "WB" uses the wild bootstrap algorithm.
}
  \item{WB_dist}{
The distribution used for the wild bootstrap. Either "rademacher", "normal", or "mammen".
}
  \item{x}{Object with class attribute \sQuote{cointBootTest}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{

\ifelse{latex}{

Consider the \eqn{K}-dimensional heteroskedastic cointegrated VAR model of Cavaliere,
Rahbek and Taylor (2014):

\deqn{
\Delta\mathbf{y}_{t}=\boldsymbol{\alpha}\boldsymbol{\beta}^{\prime}\mathbf{y}_{t-1}+\sum_{i=1}^{p-1}\mathbf{\Gamma}_{i}\Delta\mathbf{y}_{t-i}+\boldsymbol{\alpha}\boldsymbol{\rho}^{\prime}D_{t}+\bm{\phi}d_{t}+\bm{\varepsilon}_{t},\ \ \ \ t=1,\ldots ,T,
}

where \eqn{\boldsymbol{\alpha}}  and \eqn{\boldsymbol{\beta}}   are  \eqn{\left( K \times r \right)}  matrices of rank  \eqn{r<K}, the number  \eqn{r}  being the cointegration rank. \eqn{D_{t}} and \eqn{d_{t}} differ according to the \code{model} argument in the following manner:

model 1: \eqn{D_{t}=0} and \eqn{d_{t}=0} (no deterministic
terms) \cr 
model 2: \eqn{D_{t}=1} and \eqn{d_{t}=0} (restricted constant) \cr
model 3: \eqn{D_{t}=1} and \eqn{d_{t}=1} (restricted linear trend) \cr


The likelihood ratio (LR) statistic for testing cointegration rank \eqn{r} against \eqn{K} is

\deqn{
Q_{r,T}=-T\sum_{i=r+1}^{K}\mbox{log}( 1 - \widehat{\lambda}_{i},),
}

where the eigenvalues \eqn{\mathbf{\widehat{\lambda}}_{1}>\ldots >\mathbf{\widehat{\lambda}}_{K}} are the \eqn{K} largest solutions to a certain eigenvalue
problem (see Johansen 1996). \cr 


\bold{Bootstrap and wild bootstrap algorithm of Cavaliere, et al. (2012, 2014):}

\bold{1.} Estimate the model under  \eqn{H(r)} using Gaussian PMLE yielding the estimates \eqn{\boldsymbol{\widehat{\beta}}^{(r)}}, \eqn{\boldsymbol{\widehat{\alpha}}^{(r)}},
\eqn{\boldsymbol{\widehat{\rho}}^{(r)}}, \eqn{\mathbf{\widehat{\Gamma}}_{1}^{(r)},\ldots,\mathbf{\widehat{\Gamma}}_{p-1}^{(r)}} and \eqn{\widehat{\bm{\phi}}^{(r)}}, together with the corresponding residuals, \eqn{\widehat{\bm{\varepsilon}}_{r,t}}. \cr
\bold{2.} Check that the equation \eqn{| \widehat{\mathbf{A}}^{(r)}(z) |=0}, with \eqn{\widehat{\mathbf{A}}^{(r)}(z):=(1-z)\mathbf{I}_{K}-\boldsymbol{\widehat{\alpha}}^{(r)}\boldsymbol{\widehat{\beta}}^{(r)\prime}z-\sum_{i=1}^{p-1}\mathbf{\widehat{\Gamma}}_{i}^{(r)}(1-z)z^{i}}, has \eqn{K-r} roots equal to 1 and all other roots outside the unit circle. If so, procede to step 3. \cr
\bold{3.} Construct the bootstrap sample recursively from 
\deqn{
\Delta\mathbf{y}_{r,t}^{\ast}=\boldsymbol{\widehat{\alpha}}^{(r)}\boldsymbol{\widehat{\beta}}^{(r)\prime}\mathbf{y}_{r,t-1}^{\ast}+\sum_{i=1}^{p-1}\mathbf{\widehat{\Gamma}}_{i}^{(r)}\Delta\mathbf{y}_{r,t-i}^{\ast}+\boldsymbol{\widehat{\alpha}}^{(r)}\boldsymbol{\widehat{\rho}}^{(r)\prime}D_{t}+\widehat{\bm{\phi}}^{(r)}d_{t}+\bm{\varepsilon}_{r,t}^{\ast},\ \ \ \ t=1,\ldots ,T,
}
initialized at \eqn{\mathbf{y}_{r,j}^{\ast}=\mathbf{y}_{j},  j=1-p,\ldots,0}, and with the \eqn{T} bootstrap errors \eqn{\bm{\varepsilon}_{r,t}^{\ast}} generated using the residuals \eqn{\widehat{\bm{\varepsilon}}_{r,t}}. The bootstrap errors are generated depending on the \code{boot_type} argument in the following manner: \cr
\code{boot_type = "B"}: The i.i.d. bootstrap, such that \eqn{\bm{\varepsilon}_{r,t}^{\ast}:=\mathbf{\widehat{\bm{\varepsilon}}}_{r,\mathcal{U}_{t}}}, where \eqn{\mathcal{U}_{t}, t=1,\ldots ,T} is an i.i.d. sequence of discrete uniform distributions on \eqn{\{1,2,\ldots, T\}}.  \cr
\code{boot_type = "WB"}: The wild bootstrap, where for each \eqn{t=1,\ldots ,T},   \eqn{\bm{\varepsilon}_{r,t}^{\ast}:=\widehat{\bm{\varepsilon}}_{r,t}w_{t}}, where \eqn{w_{t}, t=1,\ldots ,T}, is an i.i.d. sequence distributed according to the \code{WB_dist} argument. \cr
\bold{4.} Using the bootstrap sample, \eqn{\{\mathbf{y}_{r,t}^{\ast}\}}, and denoting by \eqn{\mathbf{\widehat{\lambda}}_{1}^{\ast}>\ldots >\mathbf{\widehat{\lambda}}_{K}^{\ast}} the ordered solutions to the bootstrap analogue of the eigenvalue problem, compute the bootstrap LR statistic \eqn{Q_{r,T}^{\ast}:=-T\sum_{i=r+1}^{K}\mbox{log}( 1 - \widehat{\lambda}_{i}^{\ast})}. Define the corresponding \eqn{p}-value as \eqn{p_{r,T}^{\ast}:=1-G_{r,T}^{\ast}(Q_{r,T})}, \eqn{G_{r,T}^{\ast}}(.) denoting the conditional (on the origional data) cdf of \eqn{Q_{r,T}^{\ast}}. \cr
\bold{5.} The bootstrap test of \eqn{H(r)} against \eqn{H(K)} at level \eqn{\mathbf{\eta}} rejects \eqn{H(r)} if \eqn{p_{r,T}^{\ast}\leq\mathbf{\eta}}. \cr

If \code{r = "sequence"}, the algorithm is repeated for each null hypothesis \eqn{H(r),  r=0,\ldots ,K-1}, and the first null hypothesis with a \eqn{p_{r,T}^{\ast}>\mathbf{\eta}} is selected as the cointegration rank. If \eqn{p_{r,T}^{\ast}\leq\mathbf{\eta}, r=0,\ldots ,K-1}, the rank selected is \eqn{\widehat{r}=K}.


}{Please see the pdf version of the manual at the package's CRAN page for mathematical details of the test.}

}


\value{

a list of class \code{"cointBootTest"}. 

\item{eigen_val}{
        the eigenvalues.
}
\item{eigen_vec}{
        the eigenvectors.
}
\item{alpha}{
        a matrix with the estimated alpha parameters for the model with \code{r = K} (for other values of \code{r}, the alpha parameters are the first \code{r} columns of this matrix). 
}
\item{beta}{
        a matrix with the estimated beta parameters for the model with \code{r = K} (for other values of \code{r}, the beta parameters are the first \code{r} columns of this matrix). 
}
\item{gamma}{
        a list of matrices with the estimated gamma parameters. Each parameter matrix corresponds to the model estimated under the null hypothesis in \code{r} (0:(K-1) if \code{r = "sequence"}), in the same order. 
}
\item{rho}{
        a matrix with the estimated rho parameters for the model with \code{r = K} (for other values of \code{r}, the rho parameters are the first \code{r} columns of this matrix). 
}
\item{phi}{
        a list of matrices with the estimated phi parameters. Each parameter matrix corresponds to the model estimated under the null hypothesis in \code{r} (0:(K-1) if \code{r = "sequence"}), in the same order. 
}
\item{dummy_coefs}{
        a list of matrices with the estimated dummy parameters. Each parameter matrix corresponds to the model estimated under the null hypothesis in \code{r} (0:(K-1) if \code{r = "sequence"}), in the same order. 
}
\item{residuals}{
        a list of residual matrices, one for each model estimated under the null hypothesis in \code{r} (0:(K-1) if \code{r = "sequence"}), in that order. 
}
\item{Q}{
        a vector with the Q test statistics. If \code{r = "sequence"}, then the first element is for the null hypothesis 
        r = 0, and the last is for r = K - 1. Otherwise, the order corresponds to the  \code{r} argument.
}
\item{B.Q}{
        a matrix of the iid bootstrap Q statistics. Each column represent the null hypothesis in the order of \code{r} 
        (0:(K-1) if \code{r = "sequence"}). 
}
\item{WB.Q}{
        a matrix of the wild bootstrap Q statistics. Each column represent the null hypothesis in the order of \code{r} 
        (0:(K-1) if \code{r = "sequence"}). 
}
\item{B.r}{
        the selected cointegration rank from the iid bootstrap test, if \code{r = "sequence"}) were used.
}
\item{WB.r}{
        the selected cointegration rank from the wild bootstrap test, if \code{r = "sequence"}) were used.
}
\item{B.pv}{
        a vector with the bootstrap P.values, in the order of \code{r} (0:(K-1) if \code{r = "sequence"}). 
}
\item{WB.pv}{
        a vector with the wild bootstrap P.values, in the order of \code{r} (0:(K-1) if \code{r = "sequence"}). 
}
\item{B.errors}{
        the number of times the bootstrap simulations had to be resimulated due to errors.
}
\item{WB.errors}{
        the number of times the wild bootstrap simulations had to be resimulated due to errors.
}
\item{companion_eigen}{
        a list of matrices with the eigenvalues of the companion matrix. The inverse of the eigenvalues are the roots in step 2 of the boostrap algorithm (see the .pdf version of this help file).
}

}

\references{

Ahlgren, N. & Catani, P. (2018).
\emph{Practical Problems with Tests of Cointegration Rank with Strong Persistence and Heavy-Tailed Errors}. In Corazza, M., Durábn, M., Grané, A., Perna, C., Sibillo, M. (eds) Mathematical and Statistical Methods for Actuarial Sciences and Finance, Cham, Springer.

Cavaliere, G., Rahbek, A., & Taylor, A. M. R. (2012). 
\emph{Bootstrap determination of the co-integration rank in vector autoregressive models}, Econometrica, 80, 1721-1740.

Cavaliere, G., Rahbek, A., & Taylor, A. M. R. (2014). 
\emph{Bootstrap determination of the co-integration rank in heteroskedastic VAR models}, Econometric Reviews, 33, 606-650.

Johansen, S. (1996). 
\emph{Likelihood-based inference in cointegrated vector autoregressive models}, Oxford, Oxford University Press.

}


\examples{
\dontrun{

test <- cointBootTest(y = VodafoneCDS, r = "sequence", p = 2, model = 3, signif = 0.05, 
  dummies = NULL, B = 999, boot_type = c("B", "WB"), WB_dist = "rademacher")
test
}
}

