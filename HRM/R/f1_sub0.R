####################################################################################################################################
### Filename:    f1_sub0.R
### Description: Function for calculating the test statistic for only one whole-plot factor
###
###
###
####################################################################################################################################

#' Test for main group effect (weighted/unweighted)
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.1.none <- function(X, alpha , group, subject, data, formula, nonparametric ){

  ranked <- NULL
  varQGlobal <- NULL

  temp0 <- hrm.1w.0f(X, alpha , group,  subject, data, "A", paste(as.character(group), " weighted"), nonparametric, ranked, varQGlobal)
  temp1 <- hrm.1w.0f(X, alpha , group,  subject, data, "Au", paste(as.character(group), " unweighted"), nonparametric, ranked, varQGlobal)

  output <- list()
  output$result <- rbind(temp0, temp1)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c("none"))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- FALSE
  rownames(output$result) <- 1:dim(output$result)[1]
  class(output) <- "HRM"

  return(output)
}


#' Test for interaction of factor A and B
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the response variable
#' @param H string specifying the hypothesis
#' @param text a string, which will be printed in the output
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1w.0f <- function(X, alpha, group, subject, data, H, text, nonparametric, ranked, varQGlobal ){

  stopifnot(is.data.frame(X),is.character(subject), is.character(group),alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0


  group <- as.character(group)
  subject <- as.character(subject)

  X <- as.data.table(X)
  setnames(X, c(data, group, subject), c("data", "group", "subject"))

  a <- nlevels(X[,group])
  d <- 1
  c <- 1
  n <- table(X[,group])/d

  if(nonparametric & is.null(ranked)) {
    X[,data:= 1/(sum(n)*d)*(pseudorank(X[,data], X[, group]) - 1/2)]
  }

  X <- split(X, X[,group], drop=TRUE)
  for(i in 1:a){
    X[[i]] <- X[[i]][ order(subject), ]
    X[[i]] <- X[[i]][,data]
    X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
  }


  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }

  # creating X_bar (list with a entries)
  X_bar <- sapply(X, colMeans, na.rm=TRUE)


  if(H=="A"){
    K <- 1
    S <- diag(n)-1/sum(n)*n%*%t(n)
  } else if(H=="Au"){
    K <- 1
    S <- P(a)
  }


  # creating dual empirical covariance matrices
  K_A <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)

  ##########################
  ### U statistics
  #########################

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,K)
    }
  }

  #################################################################################################

  # f
  f_1 <- 0
  f_2 <- 0

  for(i in 1:a){
    f_1 <- f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]],nonparametric,Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]],nonparametric,Q)
    j <- i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*K%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%K%*%t(X[[i]])%*%P(n[i]))
      j <- j+1
    }
  }

  f <- f_1/f_2


  ##################################################################################################



  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0


  for(i in 1:a){
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0 <- f0_1/f0_2

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

  test <- (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  if(nonparametric) {
    output$np.correction <- FALSE
  }

  return (output)
}

# End ------------------------------------------------------------
