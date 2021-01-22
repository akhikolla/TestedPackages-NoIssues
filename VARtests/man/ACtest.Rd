\name{ACtest}

\alias{ACtest}
\alias{print.ACtest}

\title{
Test for Error Autocorrelation in VAR Models.
}
\description{
Performs LM tests for error AC in VAR models. The code is based on Paul Catani's original R code used in the paper \emph{Wild Bootstrap Tests for Autocorrelation in Vector Autoregressive Models} (Ahlgren and Catani, 2016).
}

\usage{
ACtest(fit, h = 4, HCtype = c("LM", "HC0", "HC1", "HC2", "HC3"), univariate = FALSE)
\method{print}{ACtest}(x, ...)
}



\arguments{

  \item{fit}{
an object of class \code{"VARfit"} that was returned by the \code{\link{VARfit}} function, or an object of class \code{"verest"} from the function \code{"VAR"} in the \code{"vars"} package.
}
  \item{h}{
the lag length of the alternative VAR(h) model for the errors. 
}
  \item{HCtype}{
a vector containing some or all (default) of \code{"LM"}, \code{"HC0"}, \code{"HC1"}, \code{"HC2"} and \code{"HC3"}. The \code{"LM"} type assumes homoscedastic errors, while the \code{"HC"} types are different Heteroskedasticity-consistent tests (see 'Details').
}
  \item{univariate}{
either \code{FALSE}, \code{TRUE}, or \code{"only"}. If \code{TRUE}, the test will also be performed in its univariate form, if \code{"only"}, only the univariate version will be used.
}
  \item{x}{Object with class attribute \sQuote{ACtest}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{

\ifelse{latex}{

To run the wild bootstrap version of the test, please use the output from this function with the function \code{\link{wildBoot}}.

\bold{Mathematical details}

The tests for error AC are based on the least squares (LS) residuals from the
\eqn{K}-dimensional vector autoregressive (VAR) model with \eqn{p} lags (abstracting from
deterministic terms):

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
_{1},\ldots ,\mathbf{\Pi }_{p}}.

The LM statistic is computed from the auxilary model

\deqn{
\widehat{\mathbf{u}}_{t}=\mathbf{\Pi }_{1}\mathbf{y}
_{t-1}+\cdots +\mathbf{\Pi }_{p}\mathbf{y}_{t-p} + \mathbf{D}_{1}\widehat{\mathbf{u}}_{t-1}+\cdots +%
\mathbf{D}_{h}\widehat{\mathbf{u}}_{t-h}+\mathbf{e}_{t}
}
\deqn{
=(\mathbf{Z}_{t-1}^{\prime }\otimes \mathbf{I}_{K})\bold{\phi }+(\widehat{%
\mathbf{U}}_{t-1}^{\prime }\otimes \mathbf{I}_{K})\bold{\psi }+\mathbf{e}%
_{t}, \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
}
where \eqn{\mathbf{Z}_{t-1}=(\mathbf{y}_{t-1}^{\prime },\ldots ,\mathbf{y}_{t-p}^{\prime })^{\prime }}, \eqn{\bold{\phi }
=\mbox{vec}(\mathbf{Pi}_{1},\ldots ,\mathbf{Pi}_{p})^{\prime }}, \eqn{\widehat{\mathbf{U}}%
_{t-1}=(\widehat{\mathbf{u}}_{t-1}^{\prime },\ldots ,\widehat{%
\mathbf{u}}_{t-h}^{\prime })^{\prime }} and \eqn{\bold{\psi } =\mbox{vec}(\mathbf{D}_{1},\ldots
,\mathbf{D}_{h})^{\prime }}.The symbol \eqn{\otimes } denotes the Kronecker product and
the symbol vec denotes the column vectorisation operator. The first \eqn{h}
values of the residuals \eqn{\widehat{\mathbf{u}}_{t}} are set to zero in the
auxiliary model, so that the series length is equal to the series length in
the original VAR model.

The LM statistic for error AC of order \eqn{h} is given by%
\deqn{
Q_{\mbox{LM}}(h)=N\widehat{\bold{\psi }}^{\prime }(\widehat{\boldsymbol{\Sigma 
}}^{\psi \psi })^{-1}\widehat{\bold{\psi }},  \label{LMstat}
}
where \eqn{\widehat{\bold{\psi }}} is the LS estimate of \eqn{\bold{\psi } 
} and \eqn{\widehat{\boldsymbol{\Sigma }}^{\psi \psi }} is the block of 

\deqn{
\left( N^{-1}\sum_{t=1}^{N}\left[ 
\begin{array}{c}
\mathbf{Z}_{t-1}\otimes \mathbf{I}_{K} \\ 
\widehat{\mathbf{U}}_{t-1}\otimes \mathbf{I}_{K}%
\end{array}%
\right] \widehat{\boldsymbol{\Sigma }}_{\mathbf{u}}^{-1}\left[ 
\begin{array}{cc}
\mathbf{Z}_{t-1}^{\prime }\otimes \mathbf{I}_{K} & \widehat{\mathbf{U}}%
_{t-1}^{\prime }\otimes \mathbf{I}_{K}%
\end{array}%
\right] \right) ^{-1}
}
corresponding to \eqn{\bold{\psi }}. Here \eqn{\widehat{\boldsymbol{\Sigma }}_{\mathbf{u}%
}=N^{-1}\sum_{t=1}^{N}\widehat{\mathbf{u}}_{t}\widehat{\mathbf{u}}%
_{t}^{\prime }} is the estimator of the error covariance matrix from the
VAR model.

The multivariate heteroskedasticity-consistent covariance matrix estimator (HCCME)
for the auxilary model is given by (Hafner and Herwartz 2009)
\deqn{
\mathbf{V}_{N}^{-1}\mathbf{W}_{N}\mathbf{V}_{N}^{-1}=(\boldsymbol{\Gamma }%
_{N}\otimes \mathbf{I}_{K})^{-1}\mathbf{W}_{N}(\boldsymbol{\Gamma }%
_{N}\otimes \mathbf{I}_{K})^{-1},  \label{HCCME}
}
where%
\deqn{
\mathbf{V}_{N}  =\boldsymbol{\Gamma }_{N}\otimes \mathbf{I}_{K},
}
\deqn{
\boldsymbol{\Gamma }_{N} =\frac{1}{N}\sum_{t=1}^{N}\left( 
\begin{array}{c}
\widehat{\mathbf{U}}_{t-1} \\ 
\mathbf{Z}_{t-1}%
\end{array}%
\right) \left( 
\begin{array}{cc}
\widehat{\mathbf{U}}_{t-1}^{\prime } & \mathbf{Z}_{t-1}^{\prime }%
\end{array}%
\right) ,
}
\deqn{
\mathbf{W}_{N} =\frac{1}{N}\sum_{t=1}^{N}\left( 
\begin{array}{c}
\widehat{\mathbf{U}}_{t-1} \\ 
\mathbf{Z}_{t-1}%
\end{array}%
\right) \left( 
\begin{array}{cc}
\widehat{\mathbf{U}}_{t-1}^{\prime } & \mathbf{Z}_{t-1}^{\prime }%
\end{array}%
\right) \otimes (\widehat{\mathbf{u}}_{t}\widehat{\mathbf{u}}_{t}^{\prime }).
}

The HCCME-based LM statistics for error AC are obtained from the expression for \eqn{Q_{\mbox{LM}}(h)}
by replacing \eqn{\widehat{\boldsymbol{\Sigma }}^{\psi \psi }} by the block of \eqn{%
\mathbf{V}_{N}^{-1}\mathbf{W}_{N}\mathbf{V}_{N}^{-1}=(\Gamma _{N}\otimes \mathbf{I}_{K})^{-1}\mathbf{W}_{N}(\Gamma
_{N}\otimes \mathbf{I}_{K})^{-1}} corresponding to \eqn{\bold{\psi }} and with \eqn{%
\widehat{\mathbf{u}}_{t}} defined by \eqn{HC_{0}}, \eqn{HC_{1}}, \eqn{HC_{2}} and \eqn{HC_{3}
}, respectively. \cr

\eqn{HC_{0}} uses \eqn{\widehat{\mathbf{u}}_{t}\widehat{\mathbf{u}}_{t}^{\prime }}. \cr \cr
\eqn{HC_{1}} multiplies the elements of \eqn{\widehat{\mathbf{u}}_{t}\widehat{\mathbf{u}}_{t}^{\prime }} by 
\eqn{N/(N-Kp)}. \cr \cr
\eqn{HC_{2}} replaces \eqn{\widehat{%
\mathbf{u}}_{t}} by \eqn{\widehat{\mathbf{u}}_{t}/(1-h_{t})^{1/2}}, where \eqn{h_{t}=\mathbf{Z}_{t}(\mathbf{Z}^{\prime }\mathbf{Z})^{-1}\mathbf{Z}_{t}^{\prime }} is the \eqn{t}th
diagonal element of \eqn{\mathbf{Z}(\mathbf{Z}^{\prime }\mathbf{Z})^{-1}\mathbf{Z}^{\prime }}, and \eqn{\mathbf{%
Z}=(\mathbf{Z}_{0},\ldots ,\mathbf{Z}_{N-1})}.\cr \cr
\eqn{HC_{3}} replaces \eqn{\widehat{\mathbf{u}}_{t}} by \eqn{\widehat{\mathbf{u}}_{t}/(1-h_{t})}. \cr \cr

See MacKinnon and White (1985) for details.

The recursive-design wild bootstrap (WB) tests for error AC are computed using
Algorithm 1 in Ahlgren and Catani (2016). The Fixed-design WB tests for error AC
are computed using Algorithm 2 in Ahlgren and Catani (2016).

}{To run the wild bootstrap version of the test, please use the output from this function with the function \code{\link{wildBoot}}.

For mathematical details of the test, see the pdf version of the manual at the package's CRAN web page.}


}


\value{

a list of class \code{"ACtest"}. 

\item{fit}{
        the \code{fit} argument object.
}\item{inputType}{
        the type of object of \code{fit}.
}
\item{HCtype}{
        a vector of the \code{HCtype}'s used.
}
\item{h}{
        the lag length of the alternative VAR(h) model for the errors. 
}
\item{pValues}{
        a 1 x 5 matrix of the P. values of the tests. 
}
\item{Q}{
        a 1 x 5 matrix of the Q statistics of the tests. 
}
\item{unipValues}{
        a K x 5 matrix of the P. values of the univariate tests. 
    }
\item{uniQ}{
        a K x 5 matrix of the Q statistics of the univariate tests. 
    }
\item{univariate}{
        the 'univariate' argument.
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
Ahlgren, N. & Catani, P. (2016).
\emph{Wild bootstrap tests for autocorrelation in vector autoregressive models.} Stat Papers, <doi:10.1007/s00362-016-0744-0>.

Hafner, C. M. and Herwartz, H., (2009).
\emph{Testing for Linear
Vector Autoregressive Dynamics under Multivariate Generalized Autoregressive Heteroskedasticity.} Stat Neerl, 63, 294--323

MacKinnon, J. G. and White, H. (1985).
\emph{Some Heteroskedasticity Consistent Covariance Matrix Estimators with Improved Finite Sample Properties.} J Econom, 29, 305--325

}


\seealso{
\code{\link{VARfit}} to estimate a VAR(p), and \code{\link{wildBoot}} to run the Wild Bootstrap versions of the tests.
}

\examples{

fit <- VARfit(y = VodafoneCDS, p = 3, const = TRUE, trend = FALSE)
test <- ACtest(fit = fit, h = 1, HCtype = c("LM", "HC0", "HC1", "HC2", "HC3"))

}

