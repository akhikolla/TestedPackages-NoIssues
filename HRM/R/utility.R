####################################################################################################################################
### Filename:    utility.R
### Description: trace estimator functions, dual empirical covariance matrix function
###
###
###
####################################################################################################################################

calcU <- function(X, n, grp, K){
  dat <- X[[grp]]
  if(dim(dat)[2] == 1) {
    Q <- calcUCppV(dat, n[grp], mean(dat[, 1]) )
  } else {
    Q <- calcUCpp(dat,n[grp],K,colMeans(dat))
  }
  return(Q)
}


calcU_onegroup <- function(X, n, K){
  dat <- X
  Q <- calcUCpp(dat,n,K,colMeans(X))
  return(Q)
}


#' Unbiased estimator
#'
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E1 <- function(n,i, M, nonparametric, Q) {

  trace_estimator <- ifelse(nonparametric, Q[i,1], (n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
  return (trace_estimator)
}
#' Unbiased estimator
#'
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E2 <- function(n,i, M, nonparametric, Q) {
  trace_estimator <- ifelse(nonparametric, Q[i,2], (n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
  return (trace_estimator)
}
#' Unbiased estimator
#'
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E3 <- function(M_i, M_j) {
  return (matrix.trace(M_i)*matrix.trace(M_j))
}
#' Unbiased estimator
#'
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E4 <- function(M_i,M_j) {
  return (matrix.trace(M_i%*%M_j))
}


#' Function for the output: significant p-values have on or more stars
#'
#' @param value p-value
#' @keywords internal
.hrm.sigcode <- function(value) {

  char <- ""
  if(value <= 0.1 & value > 0.05) { char <- "."}
  if(value <= 0.05 & value > 0.01) {char <- '*'}
  if(value <= 0.01 & value > 0.001) {char <- "**"}
  if(value <= 0.001 & value >= 0) {char <- "***"}
  return (char)
}


#' Function for the indentity matrix
#'
#' @param size dimension of the matrix
#' @keywords internal
I <- function(size){
  return (diag(rep(1,size)))
}
#' Function for a matrix with entries 1
#'
#' @param size dimension of the matrix
#' @keywords internal

J <- function(size){
  return (rep(1,size)%*%t(rep(1,size)))
}
#' Function for the centering matrix
#'
#' @param size dimension of the matrix
#' @keywords internal
P <- function(size){
  return (I(size)-J(size)*1/size)
}

#' Function for the dual empirical matrix
#'
#' @param Data data.frame
#' @param B not used
#' @keywords internal
DualEmpirical <- function(Data, B){
  n <- dim(Data)[1]
  B <- B(dim(Data)[2]) # B = e.g. t(J_d)%*%J_d
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}
#' Function for the dual empirical matrix
#'
#' @param Data data.frame
#' @param B part of the hypothesis matrix
#' @keywords internal
DualEmpirical2 <- function(Data, B){
  n <- dim(Data)[1]
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}
