####################################################################################################################################
### Filename:    f4_sub3.R
### Description: Function for calculating the test statistic for one whole- and three subplot factors
###
###
###
####################################################################################################################################

#' Test for 1 wholeplot and 3 subplot-factors
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable crossed with the first one
#' @param factor3 column name of the data frame X of the third factor variable crossed with the first one
#' @param data column name of the data frame X containing the observed values
#' @param S Matrix for the wholeplot factor
#' @param K1 Matrix for the first subplot factor
#' @param K2 Matrix for the second subplot factor
#' @param K3 Matrix for the third subplot factor
#' @param hypothesis String which is printed in the console as ouput to indicate which hypothesis is tested.
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1w.3f <- function(X, alpha, group , factor1, factor2, factor3, subject, data, S, K1, K2, K3, hypothesis, nonparametric, ranked, varQGlobal, np.correction ){
  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1),is.character(factor2), is.character(factor3),alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0


  group <- as.character(group)
  factor1 <- as.character(factor1)
  factor2 <- as.character(factor2)
  factor3 <- as.character(factor3)
  subject <- as.character(subject)


  X <- as.data.table(X)
  setnames(X, c(data, group, factor1, factor2, factor3, subject), c("data", "group", "factor1", "factor2", "factor3", "subject"))

  a <- nlevels(X[,group])
  d <- nlevels(X[,factor1])
  c <- nlevels(X[,factor2])
  c2 <- nlevels(X[,factor3])
  n <- table(X[,group])/(d*c*c2)


  if(nonparametric & is.null(ranked)) {
    X[,data:= 1/(sum(n)*d*c*c2)*(pseudorank(X[,data], X[, group]) - 1/2)]
  }

  X <- split(X, X[,group], drop=TRUE)
  for(i in 1:a){
    X[[i]] <- X[[i]][ order(subject, factor1, factor2, factor3), ]
    X[[i]] <- X[[i]][,data]
    X[[i]] <- matrix(X[[i]],ncol=d*c*c2,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
  }


  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))


  # creating dual empirical covariance matrices
  npcriteria <- (K1 != "P" & K2 != "P" & K3 != "P")

  kdim <- 1
  if(S == "P"){
    S <- P(a)
  } else {
    S <- 1/a*J(a)
  }
  if(K1 == "P"){
    K1 <- P(d)
    kdim <- kdim*d
  } else {
    K1 <- 1/d*J(d)
  }
  if(K2 == "P"){
    K2 <- P(c)
    kdim <- kdim*c
  } else {
    K2 <- 1/c*J(c)
  }
  if(K3 == "P"){
    K3 <- P(c2)
    kdim <- kdim*c2
  } else {
    K3 <- 1/c2*J(c2)
  }

  K <- kronecker(kronecker(K1, K2), K3)
  K_phi <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,K)
    }
  }

  eval.parent(substitute(correction <- np.correction))

  if(is.na(np.correction)) {
    eval.parent(substitute(correction <- (d*c*c2 >= max(n))))
    np.correction <- (kdim >= max(n))
  }

  if(np.correction & nonparametric) {
    if(!npcriteria) {
      for(gg in 1:a) {

          tmp <- X[[gg]]%*%K
          nr <- dim(tmp)[1]
          p <- dim(tmp)[2]
          if(nr%%2 == 1){
            nr <- nr - 1
          }
          mm <- colMeans(tmp)
          g <- rep(0,nr)
          g2 <- vector("list", length = nr)
          t2 <- matrix(rep(0,p^2), ncol = p)
          for(i in 1:nr) {
            g[i] <- t(tmp[i,] - mm) %*% (tmp[i,] - mm)
            g2[[i]] <- (tmp[i,] - mm) %*% t(tmp[i,] - mm)
            t2 <- t2 + g2[[i]]
          }

          reps <- min(150, choose(nr,nr/2))
          covs <- rep(0,reps)
          g1 <- rep(0, nr/2)
          g12 <- rep(0, nr/2)

          for(i in 1:reps) {
            grp <- sample(c(rep(1,nr/2), rep(2,nr/2)))
            g1 <- g[grp == 1]
            g12 <- g[grp == 2]
            covs[i] <- cov(g1,g12)
          }

          t4 <- rep(0, nr*(nr - 1)/2)
          k <- 1
          #t2 <- matrix(rep(0,d^2), ncol = d)
          for(i in 1:nr) {
            j <- i + 1
            while(j <= nr) {
              t4[k] <- matrix.trace(g2[[i]]%*%g2[[j]])
              k <- k + 1
              j <- j + 1
            }
          }
          # for(i in 1:(nr/2)) {
          #   t2 <- t2 + g2[[i]] + g2[[(nr/2) + i]]
          # }

          corr <- mean(covs)
          corr2 <- mean(t4) - matrix.trace((1/nr*t2)*(1/nr*t2))

          tmpQ1 <-  Q[gg,1] - corr*(n[gg]^2*1/(n[gg]^2 - n[gg]))^2
          tmpQ2 <- Q[gg,2] - corr2*(n[gg]^2*1/(n[gg]^2 - n[gg]))^2
          eval.parent(substitute(tmpQ1g <- tmpQ1))
          eval.parent(substitute(tmpQ2g <- tmpQ2))


        if(tmpQ1 > 0) {
          Q[gg,1] <- tmpQ1
        }
        if(tmpQ2 > 0) {
          Q[gg,2] <- tmpQ2
        }

      }
    }
  }

  #################################################################################################

  # f
  f_1 <- 0
  f_2 <- 0

  for(i in 1:a){
    f_1 <- f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]], nonparametric, Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j<-j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]], nonparametric, Q)
    j<-i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*K%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%K%*%t(X[[i]])%*%P(n[i]))
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
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]], nonparametric, Q)
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
  eval.parent(substitute(varQGlobal <- direct))

  test <- (t(X_bar)%*%K_phi%*%X_bar)/(t(rep(1,dim(K_phi)[1]))%*%(K_phi*direct)%*%(rep(1,dim(K_phi)[1])))
  p.value<-1-pf(test,f,f0)
  output <- data.frame(hypothesis=hypothesis,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  if(nonparametric) {
    output$np.correction <- np.correction
  }

  return (output)
}

# Hypothesis BC End ------------------------------------------------------------
