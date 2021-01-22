####################################################################################################################################
### Filename:    f2_old.R
### Description: Functions for testing when using one whole- and one subplotfactor; outdated functions only used for compatibility
###              within the function hrm.test.matrices
###
###
####################################################################################################################################


#' Test for no main treatment effect (weighted version)
#'
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A.weighted <- function(n, a, d, X, alpha, nonparametric = FALSE){

  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=2)
  stopifnot(a == length(X))
  f<-0
  f0<-0
  crit<-0
  test<-0
  if(is.data.frame(X[[1]])){
    X <- lapply(X, as.matrix)
  }


  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(as.matrix(sapply(X, colMeans, na.rm=TRUE))))

  D_a <- diag(n)-1/sum(n)*n%*%t(n)
  K_A <- kronecker(D_a,J(d))

  # creating dual empirical covariance matrices
  V <- lapply(X, DualEmpirical, B = J)


  ##########################
  ### U statistics
  #########################

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,1/d*J(d))
    }
  }


  #################################################################################################
  # f <- f_1/f_2: numerator degrees of freedom
  f_1 <- 0
  f_2 <- 0

  # computation of f_1
  for(i in 1:a){
    f_1 <- f_1 + ((1-n[i]/sum(n))^2)*.E1(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_1 <- f_1 + 2*((1-n[i]/sum(n)))*((1-n[j]/sum(n)))*.E3(V[[i]],V[[j]])
      j<-j+1
    }
  }

  # computation of f_2
  for(i in 1:a){
    f_2 <- f_2 + ((1-n[i]/sum(n)))^2*.E2(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_2 <- f_2 + 2*(n[i]*n[j])/(d^2*sum(n)^2)*.E4(d/(n[i]-1)*P(n[i])%*%X[[i]],d/(n[j]-1)*J(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%J(d)%*%t(X[[i]])%*%P(n[i]))
      j<-j+1
    }
  }

  f<-f_1/f_2

  ##################################################################################################



  #################################################################################################
  # f0 = f0_1/f0_2: denumerator degrees of freedom
  # f0_1 = f_1 numerator are the same
  f0_2 <- 0

  # computation of f0_2
  for(i in 1:a){
    f0_2 <- f0_2 + ((1-n[i]/sum(n)))^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0<-f_1/f0_2

  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)

  # constructing the  Test
  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }

  test <- (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
  p.value<-1-pf(test,f,f0)
  output <- data.frame(hypothesis="A weighted",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))

  return (output)
}

# Hypothesis 1 ------------------------------------------------------------





#' Test for no main time effect
#'
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.B <- function(n, a, d, X, alpha, nonparametric = FALSE){

  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=2, a>=1)
  stopifnot(a == length(X))
  f<-0
  f0<-0
  crit<-0
  test<-0
  X <- lapply(X, as.matrix)

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))

  # creating dual empirical covariance matrices
  V <- lapply(X, DualEmpirical, B = P)

  K_B <- kronecker(J(a),P(d))

  ##########################
  ### U statistics
  #########################

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,P(d))
    }
  }


  #################################################################################################

  # f
  f_1 <- 0
  f_2 <- 0

  for(i in 1:a){
    f_1 <- f_1 + (1/n[i])^2*.E1(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_1 <- f_1 + 2*(1/n[i])*(1/n[j])*.E3(V[[i]],V[[j]])
      j<-j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (1/n[i])^2*.E2(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_2 <- f_2 + 2*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*P(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%P(d)%*%t(X[[i]])%*%P(n[i]))
      j<-j+1
    }
  }

  f<-f_1/f_2


  ##################################################################################################



  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0


  for(i in 1:a){
    f0_2 <- f0_2 + (1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0<-f0_1/f0_2

  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)

  # Test

  direct <- 1/n[1]*var(X[[1]])
  if(a>=2){
    for(i in 2:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }


  # direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  # if(a>2){
  #   for(i in 3:a) {
  #     direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
  #   }
  # }
  test <- (t(X_bar)%*%K_B%*%X_bar)/(t(rep(1,dim(K_B)[1]))%*%(K_B*direct)%*%(rep(1,dim(K_B)[1])))
  p.value<-1-pf(test,f,f0)
  output <- data.frame(hypothesis="B",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))


  return (output)
}

# Hypothesis 2 ------------------------------------------------------------







#' Test for no interaction between treatment and time
#'
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.AB <- function(n, a, d, X, alpha, nonparametric = FALSE){

  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=2, a>=2)
  stopifnot(a == length(X))
  f<-0
  f0<-0
  crit<-0
  test<-0
  X <- lapply(X, as.matrix)

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))

  # creating dual empirical covariance matrices
  V <- lapply(X, DualEmpirical, B = P)

  K_AB <- kronecker(P(a),P(d))

  ##########################
  ### U statistics
  #########################

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,P(d))
    }
  }

  #################################################################################################
  # f
  f_1 <- 0
  f_2 <- 0


  # phi <- A
  for(i in 1:a){
    f_1 <- f_1 + (1/n[i]*(1-1/a))^2*.E1(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_1 <- f_1 + 2*(1/n[i]*(1-1/a))*(1/n[j]*(1-1/a))*.E3(V[[i]],V[[j]])
      j<-j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (1/n[i]*(1-1/a))^2*.E2(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_2 <- f_2 + 2/(n[i]*n[j]*a^2)*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*P(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%P(d)%*%t(X[[i]])%*%P(n[i]))
      j<-j+1
    }
  }

  f<-f_1/f_2


  ##################################################################################################



  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0

  for(i in 1:a){
    f0_2 <- f0_2 + (1/n[i]*(1-1/a))^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0<-f0_1/f0_2

  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)

  # Test

  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test <- (t(X_bar)%*%K_AB%*%X_bar)/(t(rep(1,dim(K_AB)[1]))%*%(K_AB*direct)%*%(rep(1,dim(K_AB)[1])))
  p.value<-1-pf(test,f,f0)
  output <- data.frame(hypothesis="AB",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))


  return (output)
}

# Hypothesis 3 ------------------------------------------------------------






#' Test for no simple treatment effect
#'
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A_B <- function(n, a, d, X, alpha, nonparametric = FALSE){

  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=2, a>=2)
  stopifnot(a == length(X))
  f<-0
  f0<-0
  crit<-0
  test<-0
  X <- lapply(X, as.matrix)

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(as.matrix(sapply(X, colMeans, na.rm=TRUE))))

  # creating dual empirical covariance matrices
  V <- lapply(X, DualEmpirical, B = I)

  K_A_B <- kronecker(P(a),I(d))

  ##########################
  ### U statistics
  #########################

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,I(d))
    }
  }

  #################################################################################################
  # f
  f_1 <- 0
  f_2 <- 0

  for(i in 1:a){
    f_1 <- f_1 + (1/n[i]*(1-1/a))^2*.E1(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_1 <- f_1 + 2*(1/n[i]*(1-1/a))*(1/n[j]*(1-1/a))*.E3(V[[i]],V[[j]])
      j<-j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (1/n[i]*(1-1/a))^2*.E2(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_2 <- f_2 + 2/(n[i]*n[j]*a^2)*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*t(X[[j]])%*%P(n[j])%*%X[[j]]%*%t(X[[i]])%*%P(n[i]))
      j<-j+1
    }
  }

  f<-f_1/f_2


  ##################################################################################################



  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0


  for(i in 1:a){
    f0_2 <- f0_2 + (1/n[i]*(1-1/a))^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0<-f0_1/f0_2

  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)

  # Test

  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test <- (t(X_bar)%*%K_A_B%*%X_bar)/(t(rep(1,dim(K_A_B)[1]))%*%(K_A_B*direct)%*%(rep(1,dim(K_A_B)[1])))
  p.value<-1-pf(test,f,f0)
  output <- data.frame(hypothesis="A|B",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))

  return (output)
}

# Hypotheses 4 ------------------------------------------------------------




#' Test for no main treatment effect (unweighted version)
#'
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A.unweighted <- function(n, a, d, X, alpha, nonparametric = FALSE){

  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=2)
  stopifnot(a == length(X))
  f<-0
  f0<-0
  crit<-0
  test<-0
  X <- lapply(X, as.matrix)

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(as.matrix(sapply(X, colMeans, na.rm=TRUE))))

  # creating dual empirical covariance matrices, where observations are transformed by a matrix B
  V <- lapply(X, DualEmpirical, B = J)

  K_A <- kronecker(P(a),1/d*J(d))

  ##########################
  ### U statistics
  #########################

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,J(d))
    }
  }


  #################################################################################################
  # f
  f_1 <- 0
  f_2 <- 0



  for(i in 1:a){
    f_1 <- f_1 + ((1-1/a)/(d*n[i]))^2*.E1(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_1 <- f_1 + 2*((1-1/a)/(d*n[i]))*((1-1/a)/(d*n[j]))*.E3(V[[i]],V[[j]])
      j<-j+1
    }
  }
  for(i in 1:a){
    f_2 <- f_2 + ((1-1/a)/(d*n[i]))^2*.E2(n,i,V[[i]],nonparametric,Q)
    j<-i+1
    while(j<=a){
      f_2 <- f_2 + 2*(1/(a^2*n[i]*n[j]*d^2))*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*J(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%J(d)%*%t(X[[i]])%*%P(n[i]))
      j<-j+1
    }
  }

  f<-f_1/f_2


  ##################################################################################################


  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0

  for(i in 1:a){
    f0_2 <- f0_2 + ((1-1/a)/(d*n[i]))^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0<-f0_1/f0_2

  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)

  # Test

  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test <- (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
  p.value<-1-pf(test,f,f0)
  output <- data.frame(hypothesis="A unweighted",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(1-pf(test,f,f0)))



  return (output)
}

# Hypothesis 1/2 ------------------------------------------------------------
