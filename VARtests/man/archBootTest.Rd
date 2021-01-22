\name{archBootTest}

\alias{archBootTest}
\alias{print.archBootTest}

\title{
Combined LM test for ARCH errors in VAR models.
}

\encoding{UTF-8}

\description{

Performs the bootstrap combined Lagrange multiplier (LM) test for autoregressive conditional heteroskedastic (ARCH) errors in vector autoregressive (VAR) models of Catani and Ahlgren (2016).

The tests of Eklund and \enc{Teräsvirta}{Terasvirta} (2007), as well as the Multivariate LM test for ARCH as described for example in \enc{Lütkepohl}{Lutkepohl} (2006, sect. 16.5), are also included if the arguments \code{ET} respectively \code{MARCH} are set to \code{TRUE}. The bootstrap procedure for those are the same as in Catani and Ahlgren (2016).

}
\usage{
archBootTest(fit, h = 2, B = 499, CA = TRUE, ET = TRUE, MARCH = TRUE,
             dist = "norm", skT.param = c(0, 1, 0, 5))
\method{print}{archBootTest}(x, ...)
}


\arguments{

  \item{fit}{
an object of class \code{"VARfit"} that was returned by the \code{\link{VARfit}} function, or an object of class \code{"verest"} from the function \code{\link[vars:VAR]{VAR}} in the \code{vars} package.
}
  \item{h}{
the lag length of the alternative VAR(h) model for the errors. 
}
  \item{B}{
the number of bootstrap simulations.
}
  \item{CA}{
if \code{TRUE}, the Catani and Ahlgren (2017) test will run.
}
  \item{ET}{
if \code{TRUE}, the Eklund and \enc{Teräsvirta}{Terasvirta} (2007) test will run.
}
  \item{MARCH}{
if \code{TRUE}, the Multivariate LM test for ARCH will run. See e.g. \enc{Lütkepohl}{Lutkepohl} (2006, sect. 16.5).
}
  \item{dist}{
the error distribution. Either \code{"norm"} for the standard normal distribution, or \code{"skT"} for the skew-t distribution. The parameters of the skew-t distribution can be set with the \code{skT.param} argument. Can also be a function that returns random draws as an (N-p) x K matrix or a vector of length (N-p) * K.
}
  \item{skT.param}{
a vector of four parameters for the skew-t distribution in case \code{"skT"} was used for the \code{dist} argument. The function \code{\link[sn:dmst]{rmst}} is used to draw the errors and the parameters are passed as \code{skT.param = c(xi, Omega, alpha, nu)}
}
  \item{x}{Object with class attribute \sQuote{archBootTest}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{

\ifelse{latex}{

All tests for ARCH are based on Cholesky-standardised least squares (LS)
residuals from the \eqn{K}-dimensional vector autoregressive (VAR) model with \eqn{p}
lags (abstracting from deterministic terms):

\deqn{
\mathbf{y}_{t}=\mathbf{\Pi }_{1}\mathbf{y}_{t-1}+\cdots +\mathbf{\Pi }_{p}
\mathbf{y}_{t-p}+\mathbf{u}_{t},\quad \mbox{E}(\mathbf{u}_{t})=\mathbf{0}
,\quad \mbox{E}(\mathbf{u}_{t}\mathbf{u}_{t}^{\prime })=\mathbf{\Omega},\ \ \ \ t=1,\ldots ,N.
}

The LS residuals are

\deqn{
\widehat{\mathbf{u}}_{t}=\mathbf{y}_{t}-\widehat{\mathbf{\Pi }}_{1}\mathbf{y}
_{t-1}-\cdots -\widehat{\mathbf{\Pi }}_{p}\mathbf{y}_{t-p},
}
where \eqn{\widehat{\mathbf{\Pi }}_{1},\ldots ,\widehat{\mathbf{\Pi }}_{p}} are
the LS estimates of the \eqn{K\times K} parameter matrices \eqn{\mathbf{\Pi }
_{1},\ldots ,\mathbf{\Pi }_{p}}. The multivariate LS residuals are \eqn{
\widehat{\mathbf{U}}=(\widehat{\mathbf{u}}_{1},\ldots ,\widehat{\mathbf{u}}
_{K})}, which is an \eqn{N\times K} matrix. The Cholesky-standardised LS 
residuals are
\deqn{
\widetilde{\mathbf{w}}_{t}=(\mathbf{S}_{\widehat{\mathbf{U}}}^{-1})^{\prime }
\widehat{\mathbf{u}}_{t},
}
where \eqn{\mathbf{S}_{\widehat{\mathbf{U}}}} is the Cholesky factor of \eqn{N^{-1}
\widehat{\mathbf{U}}^{\prime }\widehat{\mathbf{U}}}, i.e. \eqn{\mathbf{S}_{
\widehat{\mathbf{U}}}} is the (unique) upper triangular matrix such that
\deqn{
\widehat{\mathbf{\Omega }}=\mathbf{S}_{\widehat{\mathbf{U}}}^{\prime }
\mathbf{S}_{\widehat{\mathbf{U}}},\quad \widehat{\mathbf{\Omega }}
^{-1}=(N^{-1}\widehat{\mathbf{U}}^{\prime }\widehat{\mathbf{U}})^{-1}=
\mathbf{S}_{\widehat{\mathbf{U}}}^{-1}(\mathbf{S}_{\widehat{\mathbf{U}}
}^{-1})^{\prime }.
}

The LM test for ARCH of order \eqn{h} (Engle 1982) in equation \eqn{i}, \eqn{i=1,\ldots
,K}, is a test of \eqn{H_{0}:b_{1}=\cdots =b_{h}} against \eqn{H_{1}:b_{j}\neq 0}
for at least one \eqn{j\in \{1,\ldots ,h\}} in the auxiliary regression
\deqn{
\widetilde{w}_{it}^{2}=b_{0}+b_{1}\widetilde{w}_{i,t-1}^{2}+\cdots +b_{h}
\widetilde{w}_{i,t-h}^{2}+e_{it}.
}
The LM statistic has the form
\deqn{
LM_{i}=(N-p)R_{i}^{2},
}
where \eqn{R_{i}^{2}} is \eqn{R^{2}} from the auxiliary regression for equation \eqn{i}. \cr 

The combined LM statistic (Dufour et al. 2010, Catani and Ahlgren 2016) is
given by
\deqn{
\widetilde{LM}=1-\min_{1\leq i\leq K}(p(LM_{i})),
}
where \eqn{p(LM_{i})} is the \eqn{p}-value of the \eqn{LM_{i}} statistic, derived from
the asymptotic \eqn{\chi ^{2}(h)} distibution. The test is only available as a
bootstrap test. The bootstrap \eqn{p}-value is simulated using Bootstrap
Algorithm 1 of Catani and Ahlgren (2016) if the errors are normal,
\deqn{
w_{i1},\ldots ,w_{iT}\sim \mbox{N}(0,1),
}

and Bootstrap Algorithm 2 if the errors are skew-\eqn{t}  (by setting the function argument \code{dist = "skT"}),
\deqn{
w_{i1},\ldots ,w_{iT}\sim \mbox{skT}(0,1;\lambda ,v),
}

where \eqn{\lambda } is the skewness parameter and \eqn{v} is the degrees-of-freedom
parameter of the skew-\eqn{t} distribution. These parameters can be set with the \code{skT.param} argument. \cr 

The multivariate LM test for ARCH of order \eqn{h} is a generalisation of the
univariate test, and is based on the auxiliary regression

\deqn{
\mbox{vech}(\widetilde{\mathbf{u}}_{t}\widetilde{\mathbf{u}}_{t}^{\prime })=\mathbf{b}
_{0}+\mathbf{B}_{1} \mbox{vech}(\widetilde{\mathbf{u}}_{t-1}\widetilde{\mathbf{u}}_{t-1}^{\prime })
+\cdots+\mathbf{B}_{h} \mbox{vech}(\widetilde{\mathbf{u}}_{t-h}\widetilde{\mathbf{u}}_{t-h}^{\prime })
+\mathbf{e}_{t},
}

where \eqn{\mbox{vech}} is the half-vectorisation operator. The null hypothesis is \eqn{H_{0}:\mathbf{B
}_{1}=\cdots =\mathbf{B}_{h}=\mathbf{0}} against \eqn{H_{1}:\mathbf{B}_{j}\neq 
\mathbf{0\!}} for at least one \eqn{j\in \{1,\ldots ,h\}.} The
multivariate LM statistic has the form

\deqn{
MLM=\frac{1}{2}(N-p)K(K+1)-(N-p)\mbox{tr}(\widehat{\mathbf{\Omega}}_{\mbox{vech}}\widehat{\mathbf{\Omega}}^{-1}),
}

where \eqn{\widehat{\mathbf{\Omega }}_{\mbox{vech}}} is the estimator of the
error covariance matrix from the auxiliary regression and \eqn{\widehat{\mathbf{
\Omega }}} \eqn{=N^{-1}\sum_{t=1}^{N}\widetilde{\mathbf{u}}_{t}\widetilde{
\mathbf{u}}_{t}^{\prime }} is the estimator of the error covariance matrix
from the VAR model (see \enc{Lütkepohl}{Lutkepohl} 2006, sect. 16.5). The \eqn{MLM} statistic
is asymptotically distributed as \eqn{\chi ^{2}(K^{2}(K+1)^{2}h/4)}. The test is
available as an asymptotic test using the asymptotic \eqn{\chi
^{2}(K^{2}(K+1)^{2}h/4)} distribution to derive the \eqn{p}-value, and as a
bootstrap test. The bootstrap \eqn{p}-value is simulated using Bootstrap
Algorithms 1 and 2 of Catani and Ahlgren (2016). The asymptotic validity of
the bootstrap multivariate LM test has not been established. \cr


The Eklund and \enc{Teräsvirta}{Terasvirta} (2007) LM test of constant error covariance
matrix assumes the alternative hypothesis is a constant conditional
correlation autoregressive conditional heteroskedasticity (CCC-ARCH) process
of order \eqn{h}: \eqn{\mathbf{H}_{t}=\mathbf{D}_{t}\mathbf{PD}_{t}}, where \eqn{\mathbf{
D}_{t}=\mbox{diag}(h_{1t}^{1/2},\ldots ,h_{Kt}^{1/2})} is a diagonal matrix of conditional
standard deviations of the errors \eqn{\{\mathbf{u}_{t}\}} and \eqn{\mathbf{P}=(\rho
_{ij})}, \eqn{i,j=1,\ldots ,K}, is a positive definite matrix of conditional
correlations. The conditional variance \eqn{\mathbf{h}_{t}=(h_{1t},\ldots
,h_{Kt})^{\prime }} is assumed to follow a CCC-ARCH\eqn{(h)} process: 

\deqn{
\mathbf{h}_{t}=\mathbf{a}_{0}+\sum_{j=1}^{h}\mathbf{A}_{j}\boldsymbol{u}
_{t-j}^{(2)},
}
where \eqn{\mathbf{a}_{0}=(a_{01},\ldots ,a_{0K})^{\prime }} is a \eqn{K}-dimensional vector of positive constants, \eqn{\mathbf{A}_{1},\ldots ,\mathbf{A}
_{h}} are \eqn{K\times K} diagonal matrices and \eqn{\boldsymbol{u}
_{t}^{(2)}=(u_{1t}^{2},\ldots ,u_{Kt}^{2})^{\prime }}. 

The null hypothesis
is \eqn{H_{0}:\mbox{diag}(\mathbf{A}_{1})=\cdots =\mbox{diag}(\mathbf{A}_{h})=\mathbf{0}} against \eqn{H_{1}:
\mbox{diag}(\mathbf{A}_{j})\neq \mathbf{0\!}} for at least one \eqn{j\in
\{1,\ldots ,h\}}. The LM statistic has the form 
\deqn{
LM_{CCC}=(N-p)\mathbf{s}(\widehat{\boldsymbol{\theta }})^{\prime }\mathbf{I}(
\widehat{\boldsymbol{\theta }})^{-1}\mathbf{s}(\widehat{\boldsymbol{\theta }}
),
}
where \eqn{\mathbf{s}(\widehat{\boldsymbol{\theta }})} and \eqn{\mathbf{I}(\widehat{
\boldsymbol{\theta }})} are the score vector and information matrix,
respectively, estimated under the null hypothesis (see Eklund and \enc{Teräsvirta}{Terasvirta}
 2007 for details). The asymptotic distribution of the \eqn{LM_{CCC}}
statistic is \eqn{\chi ^{2}(Kh)}. The test is available as an asymptotic test using
the asymptotic \eqn{\chi ^{2}(Kh)} distribution to derive the \eqn{p}-value, and as a
bootstrap test. The bootstrap \eqn{p}-value is simulated using Bootstrap
Algorithms 1 and 2 of Catani and Ahlgren (2016). The asymptotic validity of
the bootstrap \eqn{LM_{CCC}} test has not been established.
}{See the pdf version of the manual at the package's CRAN page for mathematical details of the test.}


}

\value{

a list of class \code{"ACtest"}. 

\item{fit}{
        the \code{fit} argument object.
}\item{inputType}{
        the type of object of \code{fit}.
}
\item{h}{
        the lag length h of the alternative VAR(h) model for the errors. 
}
\item{B}{
        the number of bootstrap simulations.
}
\item{K}{
        the number of series/equations in the fitted VAR model.
}
\item{CA}{
        the \code{CA} input argument.
}
\item{ET}{
        the \code{ET} input argument.
}
\item{MARCH}{
        the \code{MARCH} input argument.
}
\item{dist}{
        the \code{dist} argument object.
}
\item{standardizedResi}{
        the Cholesky-standardized residuals.
}
\item{CA_LM}{
        the combined LM statistic of Catani and Ahlgren (2016), computed as 1 - min(P(\code{CA_LMi})).
}
\item{CA_bootPV}{
        the bootstrap P. value of the combined LM test of Catani and Ahlgren (2016).
}
\item{CA_LMi}{
        the LM statistics of Catani and Ahlgren (2016) for each time series.
}
\item{CA_LMijStar}{
        an (N-p) x K matrix of the bootstrap LM statistics for each time series (columns) and bootstrap sample (rows), 
        for the Catani and Ahlgren (2016) test.
}
\item{CA_uniBootPV}{
        a vector of length K with the univariate bootstrap P. values for each time series, for the Catani and Ahlgren 
        (2016) test.
}
\item{ET_LM}{
        the LM statistic of the Eklund and \enc{Teräsvirta}{Terasvirta} (2007) test.
}
\item{ET_PV}{
        the P.value of the Eklund and \enc{Teräsvirta}{Terasvirta} (2007) LM test statistic.
}
\item{ET_bootPV}{
        the bootstrap P.value of the Eklund and \enc{Teräsvirta}{Terasvirta} (2007) test.
}
\item{ET_LMStar}{
        the bootstrap LM test statistics for the Eklund and \enc{Teräsvirta}{Terasvirta} (2007) test.
}
\item{MARCH_LM}{
        the LM statistic of the Multivariate LM test for ARCH. See e.g. \enc{Lütkepohl}{Lutkepohl} (2006, sect. 16.5).
}
\item{MARCH_PV}{
        the P.value of the MARCH LM test statistic.
}
\item{MARCH_bootPV}{
        the bootstrap P.value of the MARCH test.
}
\item{MARCH_LMStar}{
        the bootstrap LM test statistics for the MARCH test.
}
\item{description}{
        who ran the test and when.
}
\item{time}{
        computation time taken to run the test.
}
\item{call}{
        how the function \code{ACtest()} was called.
}
}

\references{

Catani, P. and Ahlgren, N. (2016).
  \emph{Combined Lagrange multiplier test for ARCH in vector autoregressive models},
  Economics and Statistics, <doi:10.1016/j.ecosta.2016.10.006>.
  
Dufour, J.-M., Khalaf, L., and Beaulieu, M.-C. (2010). 
  \emph{Multivariate residual-based finite-sample tests for serial dependence and arch effects with applications to asset pricing models},
  Journal of Applied Econometrics, 25 (2010) 263--285.
  
Eklund, B. and \enc{Teräsvirta}{Terasvirta}, T. (2007). \emph{Testing constancy of the error covariance matrix in vector models},
  Journal of Econometrics, 140, 753-780.
    
Engle, R.F. (1982). \emph{Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflation},
  Econometrica, 50, 987-1007.
  
\enc{Lütkepohl}{Lutkepohl}, H. (2006), \emph{New Introduction to Multiple Time Series
    Analysis}, Springer, New York.

}

\seealso{
\code{\link{VARfit}} to estimate a VAR(p).
}

\examples{

fit <- VARfit(y = VodafoneCDS, p = 3, const = TRUE, trend = FALSE)
test <- archBootTest(fit = fit, h = 5, B = 199, CA = TRUE, ET = TRUE, MARCH = TRUE, dist = "norm")

}

