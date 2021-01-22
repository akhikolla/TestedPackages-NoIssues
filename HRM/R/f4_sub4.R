####################################################################################################################################
### Filename:    f4_sub4.R
### Description: Functions for calculating the test statistic for only four subplot factor
###
###
###
####################################################################################################################################

#' Test for two subplot factors
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param factor3 column name of the data frame X of the third factor variable
#' @param factor4 column name of the data frame X of the fourth factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.4.four<- function(X, alpha , factor1, factor2, factor3, factor4, subject, data, formula, testing = rep(1,2^4-1), nonparametric, np.correction ){

  # variable to store ranked observations
  ranked <- NULL
  varQGlobal <- NULL
  correction <- NULL

  # determin which hypotheses should be tested according to formula object
  hypo <- expand.grid(c(1,0),c(1,0),c(1,0),c(1,0))[1:(2^4-1),]
  hypo$sum <- rowSums(hypo)
  hypo <- hypo[order(hypo$sum, -hypo$Var1, -hypo$Var2, -hypo$Var3, -hypo$Var4),]
  x<-attributes(terms.formula(formula))$term.labels

  testing <- rep(FALSE, 2^4-1)
  for(i in 1:length(x)){
    l <- length(strsplit(x[i],":")[[1]])
    testing[which(hypo$sum == l)] <- TRUE
  }

  # create list for storing results
  temp <- list(NULL)
  for(i in 1:(2^4-1)){
    if(testing[i]) {
      temp[[i]] <-  as.data.table(hrm.0w.4s(X, alpha , factor1,  factor2, factor3, factor4, subject, data, i, "", nonparametric, ranked, varQGlobal, np.correction))
    } else {
      temp[[i]] <- NULL
    }
  }

  output <- list()
  output$result <- as.data.frame(rbindlist(temp))
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c("none"), c(factor1, factor2, factor3, factor4))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- correction
  class(output) <- "HRM"

  return (output)
}



#' Test for interaction of four subplot factors
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
hrm.0w.4s <- function(X, alpha , factor1, factor2, factor3, factor4, subject, data, H = 1, text ="", nonparametric, ranked, varQGlobal, np.correction ){

  stopifnot(is.data.frame(X),is.character(subject), is.character(factor1), is.character(factor2), alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0

  factor1 <- as.character(factor1)
  factor2 <- as.character(factor2)
  factor3 <- as.character(factor3)
  factor4 <- as.character(factor4)
  subject <- as.character(subject)

  X <- as.data.table(X)
  setnames(X, c(data, factor1, subject, factor2, factor3, factor4), c("data", "factor1", "subject", "factor2", "factor3", "factor4"))

  a <- 1
  d <- nlevels(X[,factor1])
  c <- nlevels(X[,factor2])
  c2 <- nlevels(X[,factor3])
  c3 <- nlevels(X[,factor4])
  n <- dim(X)[1]

  if(nonparametric & is.null(ranked)) {
    X[,data:= 1/(sum(n)*d*c*c2*c3)*(rank(X[,data], ties.method="average") - 1/2)]
  }

  for(i in 1:a){
    X <- X[ order(subject, factor1, factor2, factor3,factor4), ]
    X <- X[,data]
    X <- matrix(X,ncol=d*c*c2*c3,byrow=TRUE)
    n[i] <- dim(X)[1]
  }

  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }

  # creating X_bar (list with a entries)
  X_bar <- colMeans(X)  # as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))

  hypo <- expand.grid(c(1,0),c(1,0),c(1,0),c(1,0))[1:(2^4-1),]
  hypo$sum <- rowSums(hypo)
  hypo <- hypo[order(hypo$sum, -hypo$Var1, -hypo$Var2, -hypo$Var3, -hypo$Var4),][H,1:4]

  Kcomponents <- list(0)
  text <- ""
  dimension <- c(d,c,c2,c3)

  kdim <- 1
  factors <- c(factor1,factor2,factor3,factor4)
  for(i in 1:length(hypo)){

    if(hypo[i] == 1){
      Kcomponents[[i]] <- P(dimension[i])
      kdim <- kdim*dimension[i]
    } else {
      Kcomponents[[i]] <- 1/dimension[i]*J(dimension[i])
    }

    if(text == ""){
      text <- ifelse(hypo[i] == 1, paste(as.character(factors[i])), "")
    } else {
      text <- ifelse(hypo[i] == 1, paste(text, as.character(factors[i]), sep = " : "), text)
    }

  }
  K <- kronecker(kronecker(kronecker(Kcomponents[[1]], Kcomponents[[2]]), Kcomponents[[3]]), Kcomponents[[4]])

  S <- 1
  # creating dual empirical covariance matrices
  K_B <- kronecker(S, K)
  V <- list(DualEmpirical2(Data = X, B=K)) #lapply(X, DualEmpirical2, B=K)

  ##########################
  ### U statistics
  #########################
  Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU_onegroup(X,n,K)
    }
  }

  eval.parent(substitute(correction <- np.correction))

  if(is.na(np.correction)) {
    eval.parent(substitute(correction <- (d*c*c2*c3 >= max(n))))
    np.correction <- (kdim >= max(n))
  }

  if(nonparametric & np.correction) {
      for(gg in 1:a) {

        tmp <- X%*%K
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


  #################################################################################################

  # f
  f_1 <- 0
  f_2 <- 0

  for(i in 1:a){
    f_1 <- f_1 + (1*1/n[i])^2*.E1(n,i,V[[i]], nonparametric, Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(1*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (1*1/n[i])^2*.E2(n,i,V[[i]], nonparametric, Q)
    j <- i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X,1/(n[j]-1)*K%*%t(X)%*%P(n[j])%*%X%*%K%*%t(X)%*%P(n[i]))
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
    f0_2 <- f0_2 + (1*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]], nonparametric, Q)
  }

  f0 <- f0_1/f0_2

  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)

  # Test

  direct <- 1/n[1]*var(X)
  eval.parent(substitute(varQGlobal <- direct))

  test <- (t(X_bar)%*%K_B%*%X_bar)/(t(rep(1,dim(K_B)[1]))%*%(K_B*direct)%*%(rep(1,dim(K_B)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  colnames(output)[1] <- "hypothesis"
  if(nonparametric) {
    output$np.correction <- np.correction
  }

  return (output)
}

# End ------------------------------------------------------------

