####################################################################################################################################
### Filename:    f4_sub2.R
### Description: Function for calculating the test statistic for two whole- and two subplot factors
###
###
###
####################################################################################################################################

#' Test for influence of factor A
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor column name of the data frame X of within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param H string specifying the hypothesis
#' @param text a string, which will be printed in the output
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.2w.2f <- function(X, alpha, group , subgroup, factor1, factor2, subject, data, H, text = "", nonparametric, ranked, varQGlobal, np.correction ){

  stopifnot(is.data.frame(X),is.character(subject), is.character(data),is.character(group),is.character(subgroup),is.character(factor1),is.character(factor2),  alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0

  group <- as.character(group)
  subgroup <- as.character(subgroup)
  factor1 <- as.character(factor1)
  factor2 <- as.character(factor2)
  subject <- as.character(subject)
  data <- as.character((data))

  X <- as.data.table(X)
  setnames(X, c(data, group, subgroup, factor1, factor2, subject), c("data", "group", "subgroup", "factor1", "factor2", "subject"))

  ag <- nlevels(X[,group])
  asub <- nlevels(X[,subgroup])
  a <- ag*asub
  d <- nlevels(X[,factor1])
  c <- nlevels(X[,factor2])
  n <- vec(table(X[,group], X[,subgroup])/d)

  if(nonparametric & is.null(ranked)) {
    X$grouping <- as.factor(paste(X[, group], X[, subgroup], sep=""))
    X[,data:= 1/(sum(n)*d*c)*(pseudorank(X[,data], X[, grouping]) - 1/2)]
  }

  X <- dlply(as.data.frame(X), c("group", "subgroup"), .drop=TRUE)
  n <- rep(0,a)

  # X <- dlply(X, c(group, subgroup), .drop=TRUE)
  # d <- nlevels(X[[1]][,factor1])
  # n <- rep(0,a)
  # c <- nlevels(X[[1]][,factor2])

  for(i in 1:a){
    X[[i]] <- X[[i]][ order(X[[i]][,"subject"], X[[i]][,"factor1"], X[[i]][,"factor2"]), ]
    X[[i]]<- X[[i]][,"data"]
    X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
  }


  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))


  kdim <- 1
  # defining the hypothesis matrices
  if(H==1){ # A
    K <- kronecker(1/d*J(d), 1/c*J(c))
    S <- kronecker(P(ag), 1/asub*J(asub))
    text <- paste(as.character(group))
  } else if(H==2){ # A2
    K <- kronecker(1/d*J(d), 1/c*J(c))
    S <- kronecker(1/ag*J(ag), P(asub))
    text <- paste(as.character(subgroup))
  } else if(H==3){ # B
    K <- kronecker(P(d), 1/c*J(c))
    S <- kronecker(1/ag*J(ag), 1/asub*J(asub))
    text <- paste(as.character(factor1))
    kdim <- d
  } else if(H==4){ # C
    K <- kronecker(1/d*J(d), P(c))
    S <- kronecker(1/ag*J(ag), 1/asub*J(asub))
    text <- paste(as.character(factor2))
    kdim <- c
  } else if(H==5){ # AA2
    K <- kronecker(1/d*J(d), 1/c*J(c))
    S <- kronecker(P(ag), P(asub))
    text <- paste(as.character(group),":",as.character(subgroup))
  } else if(H==6){ # AB
    K <- kronecker(P(d), 1/c*J(c))
    S <- kronecker(P(ag), 1/asub*J(asub))
    text <- paste(as.character(group),":",as.character(factor1))
    kdim <- d
  } else if(H==7){ # AC
    K <- kronecker(1/d*J(d), P(c))
    S <- kronecker(P(ag), 1/asub*J(asub))
    text <- paste(as.character(group),":",as.character(factor2))
    kdim <- c
  } else if(H==8){ # A2B
    K <- kronecker(P(d), 1/c*J(c))
    S <- kronecker(1/ag*J(ag), P(asub))
    text <- paste(as.character(subgroup),":",as.character(factor1))
    kdim <- d
  } else if(H==9){ # A2C
    K <- kronecker(1/d*J(d), P(c))
    S <- kronecker(1/ag*J(ag), P(asub))
    text <- paste(as.character(subgroup),":",as.character(factor2))
    kdim <- c
  } else if(H==10){ # BC
    K <- kronecker(P(d), P(c))
    S <- kronecker(1/ag*J(ag), 1/asub*J(asub))
    text <- paste(as.character(factor1),":",as.character(factor2))
    kdim <- d*c
  } else if(H==11){ # AA2B
    K <- kronecker(P(d), 1/c*J(c))
    S <- kronecker(P(ag), P(asub))
    text <- paste(as.character(group),":",as.character(subgroup),":",as.character(factor1))
    kdim <- d
  } else if(H==12){ # AA2C
    K <- kronecker(1/d*J(d), P(c))
    S <- kronecker(P(ag), P(asub))
    text <- paste(as.character(group),":",as.character(subgroup),":",as.character(factor2))
    kdim <- c
  } else if(H==13){ # ABC
    K <- kronecker(P(d), P(c))
    S <- kronecker(P(ag), 1/asub*J(asub))
    text <- paste(as.character(group),":",as.character(factor1),":",as.character(factor2))
    kdim <- d*c
  } else if(H==14){ # A2BC
    K <- kronecker(P(d), P(c))
    S <- kronecker(1/ag*J(ag), P(asub))
    text <- paste(as.character(subgroup),":",as.character(factor1),":",as.character(factor2))
    kdim <- d*c
  } else if(H==15){ # AA2BC
    K <- kronecker(P(d), P(c))
    S <- kronecker(P(ag), P(asub))
    text <- paste(as.character(group),":",as.character(subgroup),":",as.character(factor1),":",as.character(factor2))
    kdim <- d*c
  }

  # creating dual empirical covariance matrices
  K_Hypothesis <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)

  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,K)
    }
  }

  eval.parent(substitute(correction <- np.correction))

  if(is.na(np.correction)) {
    eval.parent(substitute(correction <- (d*c >= max(n))))
    np.correction <- (kdim >= max(n))
  }

  if(np.correction & nonparametric) {
    if(H %in% c(3,4,6:15)) {
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
        for(i in 1:nr) {
          j <- i + 1
          while(j <= nr) {
            t4[k] <- matrix.trace(g2[[i]]%*%g2[[j]])
            k <- k + 1
            j <- j + 1
          }
        }

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
      j <- j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]], nonparametric, Q)
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
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]], nonparametric, Q)
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

  test <- (t(X_bar)%*%K_Hypothesis%*%X_bar)/(t(rep(1,dim(K_Hypothesis)[1]))%*%(K_Hypothesis*direct)%*%(rep(1,dim(K_Hypothesis)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  if(nonparametric) {
    output$np.correction <- np.correction
  }

  return (output)
}

# Hypothesis A End ------------------------------------------------------------
