\name{VNM-package}
\Rdversion{7.1}
\alias{VNM-package}
\docType{package}
\title{ Tools for constructing 3-objective optimal designs}
\description{The package contains several functions useful for finding 3-objective optimal designs to estimate the shape of dose-response, the ED50
	and the MED for the 4PL model. The package also can be used for the 2PL or the 3PL models.
	
	Below is the 4PL model used in this package:
	
	\bold{4-parameter logistic(4PL) model}
	\deqn{Y=\frac{\theta_1}{1+e^{\theta_2{x}+\theta_3}}+\theta_4+\varepsilon, \varepsilon \sim N(0,\sigma^{2})}
    Here \eqn{x} is the dose in log scale, \eqn{\theta_1} is the Emax, \eqn{\theta_2} is the negative of the Hill's coefficient in the Hill's model that controls the rate of change in the response, \eqn{\theta_3} is \eqn{-\theta_2} log(ED50) and \eqn{\theta_4} is the lower limit of the response when \eqn{\theta_1 >0, \theta_2 \neq  0}, and \eqn{-\infty <} ED50 \eqn{< \infty}. The 4PL model becomes the 2PL model when \eqn{\theta_1 = 1} and \eqn{\theta_4 = 0} or the 3PL model when \eqn{\theta_4 = 0}.
    
}
	
\details{
\tabular{ll}{
Package: \tab VNM\cr
Type: \tab Package\cr
Version: \tab 1.0 \cr
Date: \tab 2014-05-15\cr
License: \tab GPL-2\cr
}

The main functions are:

\itemize{
	
\bold{MOPT:} Finds and verifies that the generated design is the multiple-objective optimal design.
	
\bold{Deff:} Computes D-efficiency of the generated design for estimating the model parameters.
	
\bold{ceff1:} Computes c-efficiency of the generated design for estimating the ED50.
	
\bold{ceff2:} Computes c-efficiency of the generated design for estimating the MED.
	
\bold{S.Weight:} Computes the optimal weights for the fixed design points.
}

}
 
\author{ Seung Won Hyun, Weng Kee Wong, and Yarong Yang  }
\references{Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model. (\emph{Journal of Statistical Software}, 83, 1-19, 2018, doi: 10.18637/jss.v083.i05.)
		
	Hyun, S. W., Wong, W.K. Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. (\emph{International Journal of Biostatistics}, 11, 253-271, 2015)
}

\keyword{package}

\examples{
##Find the 3-objective optimal design for the 4PL model
result=MOPT(LB=log(.001), UB=log(1000), P=c(0.137,1.563,.00895,-1.790),lambda=c(1/3,1/3),
delta=-1)

##Generated the 3-objective optimal design
Res.D=summary(result)

##Verification plot of the generated design
plot(result)

dose=Res.D[1,]
weight=Res.D[2,]

##Check D-efficiency of the generated design
Deff(weight, P=c(0.137,1.563,.00895,-1.790), dose, LB=log(.001), UB=log(1000))

##Check c-efficiency of the generated design for estimating the ED50
#ceff1(weight, P=c(0.137,1.563,.00895,-1.790), dose, LB=log(.001), UB=log(1000))

##Check c-efficency of the generated design for estimating the MED
#ceff2(weight, P=c(0.137,1.563,.00895,-1.790), dose, LB=log(.001), UB=log(1000), delta=-1, r=30)
}
