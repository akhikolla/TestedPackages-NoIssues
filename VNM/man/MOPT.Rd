\name{MOPT}
\alias{MOPT}
\title{Search the 3-objective optimal designs for estimating model parameters, ED50, and MED.}
\description{
   Search the 3-objective optimal designs that maximizes a weighted average of efficiencies for estimating model parameters, the ED50, and the MED under the 4PL model and verifies the optimality of the generated design using the General Equivalence Theorem. This also can be used to search the 3-objective optimal designs under the 2 or 3-parameter logistic models by setting the parameter values differently.
}  
\usage{
MOPT(LB, UB, P, lambda, delta, r, grid, epsilon, epsilon_w, verbose)
}
\arguments{
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P = (p1, p2, p3, p4), where p1 is the lower limit of the response (\eqn{\theta_4}), p2 is Emax (\eqn{\theta_1}), p3 is the ED50 (exp(-\eqn{\frac{\theta_3}{\theta_2}})) and p4 is the slope at the ED50 (-\eqn{\theta_2}). For the 4-parameter logistic model, a user needs to specify all 4 nominal values in P: P = (p1, p2, p3, p4). For the 3-parameter logistic model, a user needs to specify only the 3 nominal values, Emax, the ED50, and the slope: P = (p2, p3, p4). For the 2-parameter logistic model, a user needs to specify only the 2 nominal values, the ED50 and the slope: P = (p3, p4).}
  \item{lambda}{ A numeric vector. User-selected weights for the first two objectives. lambda = c(q1, q2), where q1, q2 represent weights for estimating model parameter and estimating the ED50 respectively. They are non-negative and q1 + q2 <= 1.}
  \item{delta}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of delta units better than the minimum dose.}
  \item{r}{ Numeric. The number of iteritions to select the initial design to search the 3-objective optimal design. Default is 10 and needed to be increased (for example, r = 30 or 50) if the searched optimal design is not a true optimal. }
  \item{grid}{ Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.}
  \item{epsilon}{ Numeric. Stopping criterion for the algorithm to search the multiple-objective optimal design. Default is 0.001.}
  \item{epsilon_w}{ Numeric. Stopping criterion for the Newton Raphson method inside of the algorithm. Default is 10^-6.}
  \item{verbose}{Set to TRUE if you want more evidence of progress while data is being processed. Set to FALSE if you want your CPU cycles to be used on analysis and not printing messages. Default is FALSE.}
}
\value{An object of class OPT.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{
	Atkinson, A.C., Donev, A.N. Tobias RD Optimal Experimental Designs with SAS. (\emph{Oxford University Press, Oxford}, 2007)
	
	Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model. (\emph{Journal of Statistical Software}, 83, 1-19, 2018, doi: 10.18637/jss.v083.i05.)
	
Hyun, S. W., Wong, W.K. Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. (\emph{International Journal of Biostatistics}, 11, 253-271, 2015)
}
\examples{
	##Model parameter values for the 4PL
	par.4PL=c(0.137,1.563,.00895,-1.790)
	
	##Find the 3-objective optimal design for the 4PL model
	Res.M=MOPT(LB=log(.001), UB=log(1000), par.4PL, lambda=c(1/3,1/3), delta=-1, r=30)
	
	##Generated the 3-objective optimal design
	summary(Res.M)
	
	##Verification plot of the generated design
	plot(Res.M)
}

