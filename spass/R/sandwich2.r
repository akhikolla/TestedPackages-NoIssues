#' @title Generate Time Series with Negative Binomial Distribution and Autoregressive Correlation Structure of Order One: NB-INAR(1)
#' @description \code{rnbinom.inar1} generates one or more independent time series following the NB-INAR(1) model. The generated data has negative binomial marginal distribution and an autoregressive covariance structure.
#'
#' @param sigma   assymptotic standard deviation for Full and subpupulation
#' @param rho     correlation coefficient of the underlying autoregressive correlation structure. Must be between 0 and 1.
#' @param theta   correlation absorption coefficient if tinepoints are farther appart
#' @param k       sample size allocation factor between groups: see 'Details'.    
#' @param Time    vector of measured timepoints
#' @param dropout vector describing the percentage of dropout in every timepoint    
#' @param Model   either 1 or 2, describing if 4-regressor or 3-regressor model was used. 
#'
#' @details
#' The generated marginal negative binomial distribution with mean \code{mu} = \eqn{\mu} and \code{size} = \eqn{\eta} has density
#' \deqn{(\mu/(\mu+\eta))^x \Gamma(x + \eta)/(\Gamma(x+1)\Gamma(\eta)) (\eta/(\mu+\eta))^\eta}
#' for \eqn{0 < \mu}, \eqn{0 < \eta} and \eqn{x=0, 1, 2, ...}.
#'
#' Within the NB-INAR(1) model, the correlation between two time points \eqn{t} and \eqn{s} for \code{rho} = \eqn{\rho} is given through
#' \deqn{\rho^|t-s|}
#' for \eqn{0 \le \rho \le 1}.
#'
#' @return \code{rnbinom.inar1} returns a matrix of dimension \code{n} x \code{tp} with marginal negative binomial
#' distribution with mean \code{mu} and dispersion parameter \code{size}, and an autoregressive correlation structure
#' between time points.
#'
#' @source \code{rnbinom.inar1} computes a reparametrization of the NB-INAR(1) model by \emph{McKenzie 1986} using code contributed by Thomas Asendorf.
#'
#' @references McKenzie Ed (1986), Autoregressive Moving-Average Processes with Negative-Binomial and Geometric Marginal Distributions. \emph{Advances in Applied Probability} Vol. 18, No. 3, pp. 679-705.
#'
#' @examples
#' set.seed(8)
#' random<-rnbinom.inar1(n=1000, size=0.6, mu=2, rho=0.8, tp=6)
#' cor(random)
#'
#' #Check the marginal distribution of time point 3
#' plot(table(random[,3])/1000, xlab="Probability", ylab="Observation")
#' lines(0:26, dnbinom(0:26, mu=2, size=0.6), col="red")
#' legend("topright",legend=c("Theoretical Marginal Distribution", "Observed Distribution"), 
#' col=c("red", "black"), lty=1, lwd=c(1,2))
#'
#' @export

sandwich2 <- function(sigma, rho ,theta , k, Time, dropout, Model){
  var.F = sigma[[1]]^2
  var.S = sigma[[2]]^2
  rho = rho
  theta = theta
  tp=length(Time)
  
  remain_i = 1-dropout
  # calculate prob to remain in study to time i and j. 
  remain_ij=remain_i %*% t(remain_i) #p_{i,j}= p_i*p_j
  diag(remain_ij)=remain_i # bzw p_{i,i}=\p_i
  
  # calculate within subject correlation. As of now equal correlation in Sub and SC are presumed
      correlation = gen_cov_cor(var = 1,rho = rho,theta = theta,Time = Time,cov = FALSE)
  
  
  if (Model == 1){
    # Model 1 has 4 regresion coeffs. A better Version will be implemented soon 
    mu_j=remain_i
    t_j=Time*mu_j
    r_j=1-1/(1+k)
    r_tj=r_j*t_j
    t_j_sq=Time^2*mu_j
    r_tj_sq=r_j*t_j_sq
    
    Bread=t(matrix(c(sum(mu_j),sum(mu_j*r_j),sum(t_j),sum(r_tj),
                     sum(mu_j*r_j),sum(mu_j*r_j),sum(r_tj),sum(r_tj),
                     sum(t_j),sum(r_tj),sum(t_j_sq),sum(r_tj_sq),
                     sum(r_tj),sum(r_tj),sum(r_tj_sq),sum(r_tj_sq)
    )
    ,nrow=4))
    
    print(Bread)
    invBread=solve(Bread)
    
    # now we need to calculate our Bread for our Sandwich
    V11=0;V12=0;V13=0;V14=0
    V21=0;V22=0;V23=0;V24=0
    V31=0;V32=0;V33=0;V34=0
    V41=0;V42=0;V43=0;V44=0
    
    timesq=Time%*%t(Time)
    for(j in 1:tp){
      for(i in 1:tp){
        pre_ji=remain_ij[j,i]*correlation[j,i]
        
        V11=V11+pre_ji
        V12=V12+pre_ji*r_j
        V13=V13+pre_ji*Time[j]
        V14=V14+pre_ji*Time[j]*r_j
        
        V21=V21+pre_ji*r_j
        V22=V22+pre_ji*r_j
        V23=V23+pre_ji*r_j*Time[j]
        V24=V24+pre_ji*r_j*Time[j]
        
        V31=V31+pre_ji*Time[i]
        V32=V32+pre_ji*r_j*Time[i]
        V33=V33+pre_ji*timesq[j,i]
        V34=V34+pre_ji*r_j*timesq[j,i]
        
        V41=V41+pre_ji*r_j*Time[i]
        V42=V42+pre_ji*r_j*Time[i]
        V43=V43+pre_ji*r_j*timesq[j,i]
        V44=V44+pre_ji*r_j*timesq[j,i]
      }
    }
    Meat=t(matrix(c(V11,V12,V13,V14,
                    V21,V22,V23,V24,
                    V31,V32,V33,V34,
                    V41,V42,V43,V44),nrow=4))
    AsympVarF = invBread%*%(var.F*Meat)%*%invBread
    AsympVarS = invBread%*%(var.S*Meat)%*%invBread
    
  } else {
    
    mu_j=remain_i
    t_j=Time*mu_j
    r_j=1-1/(1+k)
    r_tj=r_j*t_j
    t_j_sq=Time^2*mu_j
    r_tj_sq=r_j*t_j_sq
    
    # Calculate asymptotic Variance of Regression coeffs in 
    Bread=t(matrix(c(sum(mu_j),sum(t_j),sum(r_tj),
                     sum(t_j),sum(t_j_sq),sum(r_tj_sq),
                     sum(r_tj),sum(r_tj_sq),sum(r_tj_sq)
    )
    ,nrow=3))
    
    invBread=solve(Bread)
    
    V11=0;V12=0;V13=0
    V21=0;V22=0;V23=0
    V31=0;V32=0;V33=0
    
    for(j in 1:tp){
      for(i in 1:tp){
        
        pre_ji=remain_ij[j,i]*correlation[j,i]
        
        V11=V11+pre_ji
        V12=V12+pre_ji * Time[j]
        V13=V13+pre_ji * Time[j]* r_j
        
        V21=V21+pre_ji * Time[i]
        V22=V22+pre_ji * Time[j] * Time[i]
        V23=V23+pre_ji * Time[j] * Time[i]*r_j
        
        V31=V31+pre_ji * Time[i] * r_j
        V32=V32+pre_ji * Time[j] * Time[i] * r_j
        V33=V33+pre_ji * Time[j] * Time[i] * r_j
      }
    }
    
    Meat=t(matrix(c(V11,V12,V13,
                    V21,V22,V23,
                    V31,V32,V33),nrow=3))
   
    AsympVarF = invBread%*%(var.F*Meat)%*%invBread
    AsympVarS = invBread%*%(var.S*Meat)%*%invBread
  }
  return(list(AsympVarF=AsympVarF,AsympVarS=AsympVarS))
}