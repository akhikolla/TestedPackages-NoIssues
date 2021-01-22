\name{Deff}
\alias{Deff}
\title{Checking D-efficiency for estimating model parameters. }
\description{
  Obtaining D-efficiency for estimating model parameters for a given design under the 4-parameter logistic model. This also can be used to compute the D-efficiency under the 2 or 3-parameter logistic models by setting the parameter values differently.
}
\usage{
Deff(weight,P,dose,LB,UB,r,grid,epsilon,epsilon_w)
}
\arguments{
  \item{weight}{ A numeric vector. Weights for a given design. The weights represent the proportional allocations of subjects to the dose levels in a given design.}
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P = (p1, p2, p3, p4), where p1 is the lower limit of the response (\eqn{\theta_4}), p2 is Emax (\eqn{\theta_1}), p3 is the ED50 (exp(-\eqn{\frac{\theta_3}{\theta_2}})) and p4 is the slope at the ED50 (-\eqn{\theta_2}). For the 4-parameter logistic model, a user needs to specify all 4 nominal values in P: P = (p1, p2, p3, p4). For the 3-parameter logistic model, a user needs to specify only the 3 nominal values, Emax, the ED50, and the slope: P = (p2, p3, p4). For the 2-parameter logistic model, a user needs to specify only the 2 nominal values, the ED50 and the slope: P = (p3, p4).  }
  \item{dose}{ A vector. Dose levels for a given design.}
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
  \item{r}{Numeric. The number of iterations to select the initial design to search D-optimal design for estimating model parameters. Default is 10 and needed to be increased (for example, r = 30 or 50) if the searched D-optimal design is not a true optimal.}
  \item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.}
  \item{epsilon}{Numeric. Stopping criterion for the algorithm to search D-optimal design. Default is 0.001.}
  \item{epsilon_w}{Numeric. Stopping criterion for the Newton Raphson method inside of the algorithm. Default is 10^-6.}
}
\value{An object of class OPT.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model. (\emph{Journal of Statistical Software}, 83, 1-19, 2018, doi: 10.18637/jss.v083.i05.)
	
Hyun, S. W., Wong, W.K. Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. (\emph{International Journal of Biostatistics}, 11, 253-271, 2015)
}
\examples{
	##The given design
	dose=c(-6.91,-4.89,-4.18,6.91)
	weight=c(.344,.323,.162,.171)
	
	##Model parameter values for the 4PL
	par.4PL=c(0.137,1.563,.00895,-1.790)
	
	##Check D-efficiency of the given design and its verification plot
	Res.D=Deff(weight, P=par.4PL, dose, LB=log(.001), UB=log(1000))
	summary(Res.D)
	plot(Res.D)
}

