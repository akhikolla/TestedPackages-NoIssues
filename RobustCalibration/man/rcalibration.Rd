\name{rcalibration}
\alias{rcalibration}
%\alias{show.rgasp}
\alias{rcalibration-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Setting up the robust Calibration model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Setting up the Calibration model for estimating the parameters via MCMC with or without a discrepancy function. 
%The range and noise-variance ratio parameters are given and/or have been estimated.
}
\usage{
  rcalibration(design, observations, p_theta=NULL,X=matrix(0,dim(design)[1],1),
  have_trend=FALSE,simul_type=0, input_simul=NULL, output_simul=NULL,simul_nug=FALSE,
  math_model=NULL,theta_range=NULL, sd_proposal=c(rep(0.05,p_theta),
  rep(0.25,dim(design)[2]),0.25),S=10000,S_0=1000, discrepancy_type='S-GaSP',
  kernel_type='matern_5_2', tilde_lambda=1/2,  a=1/2-dim(design)[2],b=1,
  alpha=rep(1.9,dim(design)[2]), output_weights=rep(1,dim(design)[1]))
 % \S4method{show}{rgasp}(object)    
      
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{a matrix of observed inputs.
%%     ~~Describe \code{design} here~~
}
  \item{observations}{ a vector of experimental data.
%%     ~~Describe \code{response} here~~
}
  \item{p_theta}{an integer about the number of parameters, which should be specified by the user. 
%%     ~~Describe \code{response} here~~
}

  \item{X}{a matrix of the mean/trend discrepancy between the reality and math model. 
%%     ~~Describe \code{trend} here~~
}
  \item{have_trend}{a bool value meaning whether we assume a mean/trend discrepancy  function.
%%     ~~Describe \code{trend} here~~
}

  \item{simul_type}{an integer about the math model/simulator. If the simul_type is 0, it means we use RobustGaSP R package to build an emulator for emulation. If the simul_type is 1, it means the function of the math model is given by the user. When simul_type is 2 or 3, the mathematical model is the geophyiscal model for Kilauea Volcano.  If the simul_type is 2, it means it is for the ascending mode InSAR data; if the simul_type is 3, it means it is for the descending mode InSAR data.
}

  \item{input_simul}{an D x (p_x+p_theta) matrix of design for emulating the math model. It is only useful if simul_type is 0, meaning that we emulate the output of the math model.
}

  \item{output_simul}{a D dimensional vector of the math model runs on the design (input_simul). It is only useful if simul_type is 0, meaning that we emulate the output of the math model.
}
  \item{simul_nug}{a bool value meaning whether we have a nugget for emulating the math model/simulator. If the math model is stochastic, we often need a nugget. If simul_Nug is TRUE, it means we have a nugget for the emulator. If simul_Nug is FALSE, it means we do not have a nugget for the emulator.
%%     ~~Describe \code{trend} here~~
}


  \item{math_model}{a function of the math model provided by the user. It is only useful if simul_type is 1, meaning that we know the math model and it can be computed fast. If the evaluation the math model is computationally slow, one should set simul_type to be 0 to emulate the math model. 
}

  \item{theta_range}{a p_theta x 2 matrix of the range of the calibration parameters. The first column is the lower bound and the second column  is the upper bound. It should be specified by the user if the simul_type is 0. 
}


  \item{sd_proposal}{a vector of the standard deviation of the proposal distribution in MCMC.
}


  \item{S}{an integer about about how many posterior samples to run. 
}
  \item{S_0}{an integer about about the number of burn-in samples. 
}

  \item{discrepancy_type}{characters about the type of the discrepancy.  If it is 'no-discrepancy', it means no discrepancy function. If it is 'GaSP', it means the GaSP model for the discrepancy function. If it is 'S-GaSP', it means the S-GaSP model for the discrepancy function.}


  \item{kernel_type}{characters about the type of the discrepancy.type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern kernel} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential kernel with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha.}
 \item{tilde_lambda}{a numeric value about how close the math model to the reality in squared distance when the S-GaSP model is used for modeling the discrepancy.}
 \item{a}{a scalar of the prior parameter.}
  \item{b}{a scalar of the prior parameter.}
  \item{alpha}{a numeric parameter for the roughness in the kernel.}
  \item{output_weights}{a vector of the weights of the outputs.}

%  \item{post_sample}{a matrix of the posterior samples after burn-in.}
%  \item{post_value}{a vector of the posterior values after burn-in.}
%  \item{accept_S}{a vector of  the number of proposed samples of the calibation parameters are accepted in MCMC. The first value is the number of proposed calibration parameters  are accepted in MCMC. The second value is  the number of proposed range and nugget parameters  are accepted, if  \code{discrepancy_type}  is specified as 'GaSP' or 'S-GaSP'.}
%  \item{count_boundary}{a vector of  the number of proposed samples of the calibation parameters are outside the range and they are rejected directly.}


}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%expand here the details.
%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{rcalibration} returns an S4 object of class \code{rcalibration} (see \code{rcalibration-class}).

%If there is an emulator (i.e. simul_type is 0), \code{rcalibration} returns a list with 
%\item{rcalibration}{an S4 object of class \code{rcalibration} (see \code{rcalibration-class}.}
%\item{emulator}{an S4 object of class \code{rgasp} produced by RobustGaSP R Package .}


}

\references{
%% ~put references to the literature/web site here ~

A. O'Hagan and M. C. Kennedy (2001), \emph{Bayesian calibration of computer models}, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology}, \bold{63}, 425-464.

K. R. Anderson and M. P. Poland (2016), \emph{Bayesian estimation of magma supply, storage, and eroption rates using a multiphysical volcano model: Kilauea volcano, 2000-2012.}. \emph{Eath and Planetary Science Letters}, \bold{447}, 161-171.

K. R. Anderson and M. P. Poland (2017), \emph{Abundant carbon in the mantle beneath Hawaii}. \emph{Nature Geoscience}, \bold{10}, 704-708.


M. Gu (2016), \emph{Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output}, Ph.D. thesis., Duke University.

M. Gu and L. Wang (2017) \emph{Scaled Gaussian Stochastic Process for Computer Model Calibration and Prediction}. arXiv preprint arXiv:1707.08215.

M. Gu (2018) \emph{Jointly Robust Prior for Gaussian Stochastic Process in Emulation, Calibration and Variable Selection
}. arXiv preprint arXiv:1804.09329.
}

\author{
\packageAuthor{RobustCalibration}

Maintainer: \packageMaintainer{RobustCalibration}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
library(RobustCalibration)

#-------------------------------------------------------------
# an example with multiple local maximum of minimum in L2 loss
#-------------------------------------------------------------

## the reality 
test_funct_eg1<-function(x){
  x*cos(3/2*x)+x
}



## obtain 25 data from the reality plus a noise
set.seed(1)
## 10 data points are very small, one may want to add more data
n=15
input=seq(0,5,5/(n-1))
input=as.matrix(input)

output=test_funct_eg1(input)+rnorm(length(input),mean=0,sd=0.1)

num_obs=n=length(output)



## the math model 
math_model_eg1<-function(x,theta){
  sin(theta*x)+x  
}

##fit the S-GaSP model for the discrepancy
##one can choose the discrepancy_type to GaSP, S-GaSP or no discrepancy
##p_theta is the number of parameters to calibrate and user needs to specifiy 
##one may also want to change the number of posterior samples by change S and S_0
p_theta=1
model_sgasp=rcalibration(design=input, observations=output, p_theta=p_theta,simul_type=1,
                         math_model=math_model_eg1,theta_range=matrix(c(0,3),1,2)
                         ,S=10000,S_0=2000,discrepancy_type='S-GaSP')

##if the acceptance rate is too low or two high, one can adjust sd_proposal, e.g.
#model_sgasp=rcalibration(design=input, observations=output, p_theta=1,simul_type=1,
#                         sd_proposal=c(rep(0.02,p_theta),rep(0.2,dim(input)[2]),0.2)
#                         math_model=math_model_eg1,theta_range=matrix(c(0,3),1,2)
#                         ,S=10000,S_0=2000,discrepancy_type='S-GaSP')

##posterior samples of calibration parameter and value
plot(model_sgasp@post_sample[,1],type='l',xlab='num',ylab=expression(theta))   
plot(model_sgasp@post_value,type='l',xlab='num',ylab='posterior value')   


show(model_sgasp)




##one may want to fit a a model with an estimated baseline mean discrepancy by setting 
##X=matrix(1,dim(input_stack)[1],1),have_trend=TRUE

model_sgasp_with_mean=rcalibration(design=input, observations=output, p_theta=1,simul_type=1,
                                   X=matrix(1,dim(input)[1],1),have_trend=TRUE,
                                   math_model=math_model_eg1,theta_range=matrix(c(0,3),1,2),
                                   S=10000,S_0=2000,discrepancy_type='S-GaSP')

show(model_sgasp_with_mean)

##posterior samples of calibration parameter and value
plot(model_sgasp_with_mean@post_sample[,1],type='l',xlab='num',ylab=expression(theta))   
plot(model_sgasp_with_mean@post_value,type='l',xlab='num',ylab='posterior value')   


\dontrun{
  #-------------------------------------------------------------
  # an example with multiple local maximum of minimum in L2 loss
  # for combing the emulator
  #-------------------------------------------------------------
  
  ## the reality 
  test_funct_eg1<-function(x){
    x*cos(3/2*x)+x
  }
  
  ## obtain 20 data from the reality plus a noise
  set.seed(1)
  n=20
  input=seq(0,5,5/(n-1))
  input=as.matrix(input)
  
  output=test_funct_eg1(input)+rnorm(length(input),mean=0,sd=0.05)
  
  num_obs=n=length(output)
  
  ## the math model 
  math_model_eg1<-function(x,theta){
    sin(theta*x)+x  
  }
  
  ##let's build an emulator for the case if the math model is too slow
  
  # let's say we can only run the math model n_design times
  n_design=80
  
  design_simul=matrix(runif(n_design*2),n_design,2)
  design_simul[,1]=5*design_simul[,1]   ##the first one is the observed input x
  design_simul[,2]=3*design_simul[,2]   ##the second one is the calibration parameter 
  
  output_simul=math_model_eg1(design_simul[,1],design_simul[,2])
  
  
  
  ##this is a little slow compared with the previous model
  model_sgasp_emulator=rcalibration(design=input, observations=output, p_theta=1,simul_type=0, 
                                    input_simul=design_simul, output_simul=output_simul,
                                    theta_range=matrix(c(0,3),1,2),
                                    S=10000,S_0=2000,discrepancy_type='S-GaSP')
  
  ##now the output is a list
  show(model_sgasp_emulator)

  ##here is the plot
  plot(model_sgasp_emulator@post_sample[,1],type='l',xlab='num',ylab=expression(theta))   
  plot(model_sgasp_emulator@post_value,type='l',xlab='num',ylab='posterior value')   
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
