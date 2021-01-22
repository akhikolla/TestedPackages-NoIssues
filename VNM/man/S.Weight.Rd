\name{S.Weight}
\alias{S.Weight}
\title{Identify optimal weights for given dose levels.}
\description{
  Find optimal proportions of subjectives at the given dose levels to estimate the model parameters, the ED50 and the MED simultaneously under the 4-parameter logistiv model. This also can be used to find the optimal weights at the given dose under the 2 or 3-parameter logistic models by setting the parameter values differently.
}
\usage{
S.Weight(X,P,lambda,delta,epsilon_w)
}
\arguments{
  \item{X}{ A numeric vector. Given dose levels to search the optimal weights.}
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P=(p1, p2, p3, p4), where p1 is the lower limit of the response (\eqn{\theta_4}), p2 is Emax (\eqn{\theta_1}), p3 is the ED50 (exp(-\eqn{\frac{\theta_3}{\theta_2}})) and p4 is the slope at the ED50 (-\eqn{\theta_2}). For the 4-parameter logistiv model, a user needs to specify all 4 nominal values in P: P = (p1, p2, p3, p4). For the 3-parameter logistic model, a user needs to specify only the 3 nominal values, Emax, the ED50, and the slope: P = (p2, p3, p4). For the 2-parameter logistic model, a user needs to specify only the 2 nominal values, the ED50 and the slope: P = (p3, p4).}
  \item{lambda}{ A numeric vector. User-selected weights for the first two objectives. lambda = c(q1, q2), where q1,q2 represent weights for estimating model parameter and estimating the ED50 respectively. They are non-negative and q1 + q2 <= 1.}
  \item{delta}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of dt units better than the minimum dose.}
  \item{epsilon_w}{ Numeric. Stopping criterion for the Newton Raphson method to search the optimal weights for the given dose levels. Default is 10^-6.}
} 
\value{An object of class SW.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model. (\emph{Journal of Statistical Software}, 83, 1-19, 2018, doi: 10.18637/jss.v083.i05.)
	
Hyun, S. W., Wong, W.K. Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. (\emph{International Journal of Biostatistics}, 11, 253-271, 2015)
}
\examples{
	##The given dose levels
	dose=c(-6.91,-4.89,-4.18,6.91)
	
	##Model parameter values for the 4PL
	par.4PL=c(0.137,1.563,.00895,-1.790)
	
	##Find the optimal weights for the given dose levels
	Res.W=S.Weight(dose, par.4PL, lambda=c(1/3,1/3), delta=-1)
	
	##Print the obtained optimal weights, and its verification
    summary(Res.W)
}

