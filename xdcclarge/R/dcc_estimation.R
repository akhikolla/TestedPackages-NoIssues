#' This function estimates the parameters(alpha,beta) and time-varying correlation matrices(Rt) of DCC-GARCH model.
#' @param ini.para initial DCC-GARCH parameters(alpha,beta) of optimization
#' @param ht matrix of conditional variance vectors
#' @param residuals matrix of residual(de-mean) returns
#' @param method shrinkage method of unconditional correlation matrix(Cov:sample,LS:Linear Shrinkage,NLS:NonLinear Shrinkage)
#' @param ts ts how many time series are you taking(dufalut:1 latest value)
#'
#' @return time-varying correlations(Rt) and the result of estimation
#' @note Rt are vectorized values of the conditional correlation matrix(Rt) until time t(ts) for each row.
#' @export
#' @examples
#'   library(rugarch)
#'   library(xdcclarge)
#'   #load data
#'   data(us_stocks)
#'   n<-3
#'   Rtn<-log(us_stocks[-1,1:n]/us_stocks[-nrow(us_stocks),1:n])
#'   
#'   # Step 1:GARCH Parameter Estimation with rugarch
#'   spec = ugarchspec()
#'   mspec = multispec( replicate(spec, n = n) )
#'   fitlist = multifit(multispec = mspec, data = Rtn)
#'   ht<-sigma(fitlist)^2
#'   residuals<-residuals(fitlist)
#'   
#'   # Step 2:DCC-GARCH Parameter Estimation with xdcclarge
#'   DCC<-dcc_estimation(ini.para=c(0.05,0.93) ,ht ,residuals)
#'   #Time varying correlation matrix Rt at time t
#'   (Rt<-matrix(DCC$dcc_Rt,n,n))
#'   
#'   \dontrun{
#'   #If you want Rt at time t-s,then
#'   s<-10
#'   DCC<-dcc_estimation(ini.para=c(0.05,0.93) ,ht ,residuals,ts = s)
#'   matrix(DCC$cdcc_Rt[s,],n,n)
#'   }

dcc_estimation <- function(ini.para=c(0.05,0.93) ,ht ,residuals, method=c("COV","LS","NLS"), ts = 1){

  #de-garch residual returns
  stdresids <- residuals/sqrt(ht)

  flag <- match.arg(method)

  if(flag == "NLS"){
    print("non-linear shrinkage Step Now...")
    uncR <- nlshrink::nlshrink_cov(stdresids)
  }
  else if(flag == "LS"){
    print("linear shrinkage Step Now...")
    uncR <- nlshrink::linshrink_cov(stdresids)
  }
  else{
    uncR <- stats::cov(stdresids)
  }

  print("Optimization Step Now...")
  result <- dcc_optim(param=ini.para,ht=ht, residuals=residuals, stdresids=stdresids, uncR=uncR)

  print("Construction Rt Step Now...")
  dcc_Rt <- dcc_correlations(result$par, stdresids, uncR, ts)

  list(result=result,dcc_Rt=dcc_Rt)
}

