####################################################################################################################################
### Filename:    MVHRM_1w_1s.R
### Description: Function for calculating the MVHRM test statistic
###
###
###
####################################################################################################################################



#' MVHRM test
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
hrm.mv.1w.1f <- function(X, alpha, group , factor1, subject, data, variable, formula, nonparametric ){

  ntest <- 3

  # case 1w.1s.
  if(!is.null(group) & !is.null(factor1)) {
    group <- as.character(group)
    factor1 <- as.character(factor1)
    subject <- as.character(subject)
    variable <- as.character(variable)

    X <- as.data.table(X)
    setnames(X, c(data, group, factor1, subject, variable), c("data", "group", "factor1", "subject", "variable"))

    a <- nlevels(X[,group])
    t <- nlevels(X[,factor1])
    p <- nlevels(X[,variable])
    d <- p*t
    n <- table(X[,group])/d

    if(t <= 1) {
      stop("At least t = 2 repeated measurements needed.")
    }
    if (p <= 1) {
      stop("At least p = 2 variables needed.")
    }

    C1 <- J(a)/a
    C2 <- P(a)/(a-1)
    D1 <- t(as.matrix(rep(1,t)))*1/t
    D2 <- cbind(I(t-1), rep(-1,t-1))


    result <- data.frame(method = 1:ntest, fH = 1:ntest, fG = 1:ntest, test = 1:ntest, pvalue = 1:ntest)
    result <- hrm.mv.internal(X, "group" , "factor1", "subject", "data", "variable", C2, D1, nonparametric )
    result <- rbind(result, hrm.mv.internal(X, "group" , "factor1", "subject", "data", "variable", C1, D2, nonparametric ))
    result <- rbind(result, hrm.mv.internal(X, "group" , "factor1", "subject", "data", "variable", C2, D2, nonparametric ))

    output2 <- data.frame(hypothesis = c(rep(as.character(group), ntest), rep(as.character(factor1), ntest), rep(paste(as.character(group), ":", as.character(factor1)),ntest)), result = result)
    colnames(output2) <- c("hypothesis", "method", "fH", "fG", "test", "pvalue", "sign.code")
    for(r in 3:6) {
      output2[, r] <- as.double(output2[, r])
    }

    output <- list()
    output$result <- output2
    output$formula <- formula
    output$alpha <- alpha
    output$subject <- subject
    output$variable <- variable
    output$factors <- list(group, factor1)
    output$data <- X
    output$nonparametric <- nonparametric
    output$np.correction <- FALSE
    class(output) <- "HRM"

    return(output)
  } # end group != NULL

  # case 0w.1s.
  if(is.null(group) & !is.null(factor1)) {
    group <- "group"
    factor1 <- as.character(factor1)
    subject <- as.character(subject)
    variable <- as.character(variable)

    X$group <- as.factor(1)

    X <- as.data.table(X)
    setnames(X, c(data, group, factor1, subject, variable), c("data", "group", "factor1", "subject", "variable"))

    a <- nlevels(X[,group])
    t <- nlevels(X[,factor1])
    p <- nlevels(X[,variable])
    d <- p*t
    n <- table(X[,group])/d

    if(t <= 1) {
      stop("At least t = 2 repeated measurements needed.")
    }
    if (p <= 1) {
      stop("At least p = 2 variables needed.")
    }

    C1 <- J(a)/a
    #C2 <- P(a)/(a-1)
    D1 <- t(as.matrix(rep(1,t)))*1/t
    D2 <- cbind(I(t-1), rep(-1,t-1))



    result <- data.frame(method = 1:ntest, fH = 1:ntest, fG = 1:ntest, test = 1:ntest, pvalue = 1:ntest)
    result <- hrm.mv.internal(X, "group" , "factor1", "subject", "data", "variable", C1, D2, nonparametric )

    output2 <- data.frame(hypothesis = rep(as.character(factor1), ntest), result = result)
    colnames(output2) <- c("hypothesis", "method", "fH", "fG", "test", "pvalue", "sign.code")

    output <- list()
    output$result <- output2
    output$formula <- formula
    output$alpha <- alpha
    output$subject <- subject
    output$variable <- variable
    output$factors <- list(NULL, factor1)
    output$data <- X
    output$nonparametric <- nonparametric
    output$np.correction <- FALSE
    class(output) <- "HRM"

    return(output)
  } # end group == NULL

  # case 1w.0s.
  if(!is.null(group) & is.null(factor1)) {
    group <- as.character(group)
    factor1 <- "factor1"
    subject <- as.character(subject)
    variable <- as.character(variable)

    X$factor1 <- as.factor(1)

    X <- as.data.table(X)
    setnames(X, c(data, group, factor1, subject, variable), c("data", "group", "factor1", "subject", "variable"))

    a <- nlevels(X[,group])
    t <- nlevels(X[,factor1])
    p <- nlevels(X[,variable])
    d <- p*t
    n <- table(X[,group])/d

    if(a <= 1) {
      stop("At least a = 2 groups are needed.")
    }
    if (p <= 1) {
      stop("At least p = 2 variables needed.")
    }

    C1 <- J(a)/a
    C2 <- P(a)/(a-1)
    D1 <- t(as.matrix(rep(1,t)))*1/t
    #D2 <- cbind(I(t-1), rep(-1,t-1))



    result <- data.frame(method = 1:ntest, fH = 1:ntest, fG = 1:ntest, test = 1:ntest, pvalue = 1:ntest)
    result <- hrm.mv.internal(X, "group" , "factor1", "subject", "data", "variable", C2, D1, nonparametric )

    output2 <- data.frame(hypothesis = rep(as.character(group), ntest), result = result)
    colnames(output2) <- c("hypothesis", "method", "fH", "fG", "test", "pvalue", "sign.code")

    output <- list()
    output$result <- output2
    output$formula <- formula
    output$alpha <- alpha
    output$subject <- subject
    output$variable <- variable
    output$factors <- list(group, NULL)
    output$data <- X
    output$nonparametric <- nonparametric
    output$np.correction <- FALSE
    class(output) <- "HRM"

    return(output)
  } # end factor1 == NULL
}




#' MVHRM test
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
hrm.mv.internal <- function(X, group , factor1, subject, data, variable, C, D, nonparametric ){

    group <- as.character(group)
    factor1 <- as.character(factor1)
    subject <- as.character(subject)
    variable <- as.character(variable)

    X <- as.data.table(X)
    setnames(X, c(data, group, factor1, subject, variable), c("data", "group", "factor1", "subject", "variable"))

    a <- nlevels(X[,group])
    t <- nlevels(X[,factor1])
    p <- nlevels(X[,variable])
    d <- p*t
    n <- table(X[,group])/d

    X <- X[order(group, subject, variable, factor1),]

    if(nonparametric){
      X$data <- unlist(by(data.frame(data = X[, data], group = X[, group]), X[, variable], FUN = function(df){ pseudorank(df[, 1], df[, 2])}))
    }

    X <- X[order(group, subject, factor1, variable)]


    # EEG$v2 <- unlist(by(EEG[,c("value", "group")], EEG[, "variable"], FUN = function(df){ pseudorank(df[, 1], df[, 2])} ))
    # EEG$v2[1:10]
    # df <- subset(EEG, EEG$variable == 1)
    # dfr <- pseudorank(value~group,data=df)
    # dfr[1:10]

    R = data.frame(1:(p*t))
    for(i in 1:sum(n)){
      R[,i] <- X[((i-1)*d+1):(i*d) , data]
    }
    R  = as.matrix(R) # in columns are the Y_ik
    Z = kronecker(D, I(p))%*%R
    #Z = R
    Z = as.matrix(Z)

    Si <- list(0)
    Zi <- list(0)
    S <- 0
    Zbar <- 0
    m <- c(0, n)

    # calculate Zbar and covariance matrices
    for(i in 1:a) {
      Si[[i]] <- var(t(Z[,(sum(m[1:i])+1):(sum(m[1:(i+1)])) ]))
      Zi[[i]] <- as.matrix(1/n[i]*rowSums(Z[,(sum(m[1:i])+1):(sum(m[1:(i+1)])) ]))
      Zbar <- as.matrix(cbind(Zbar, Zi[[i]]))
      S <- S + 1/a*1/n[[i]]*Si[[i]]
    }
    Zbar <- Zbar[,-1]

    # calculate U-statistics for trace of covariance estimators
    Q <- data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
    for(i in 1:a) {
      Q[i,] <- calcUCpp(t(Z[,(sum(m[1:i])+1):(sum(m[1:(i+1)])) ]),n[i],I(dim(Si[[i]])[1]),colMeans(t(Z[,(sum(m[1:i])+1):(sum(m[1:(i+1)])) ])))
    }

    # calculate degrees of freedom
    fHden <- 0
    for(i in 1:a) {
      j <- i + 1
      fHden <- fHden + 1/n[i]^2*C[i,i]^2*(Q[i, 2] + Q[i, 1])
      while(j <= a){
        fHden <- fHden + 2*1/n[i]*1/n[j]*C[i,j]*C[j,i]*(matrix.trace(Si[[i]]%*%Si[[j]])+matrix.trace(Si[[i]])*matrix.trace(Si[[j]]))
        j <- j + 1
      }
    }

    fGden <- 0
    for(i in 1:a) {
      fGden <- fGden + 1/n[i]^2*1/(n[i]-1)^2*(n[i]*(1-1/n[i])^2 + n[i]*(n[i]-1)*1/n[i]^2)*(Q[i, 2] + Q[i, 1])
    }
    fGden <- fGden/a^2

    num <- 0
    for(i in 1:a) {
      j <- i + 1
      num <- num + 1/n[i]^2*(Q[i, 2] + Q[i, 1])
      while(j <= a){
        num <- num + 2*1/n[i]*1/n[j]*(matrix.trace(Si[[i]]%*%Si[[j]])+matrix.trace(Si[[i]])*matrix.trace(Si[[j]]))
        j <- j + 1
      }
    }
    num <- num/a^2

    fH <- num/fHden
    fG <- num/fGden

    # Test Statistic
    H = fH*Zbar%*%C%*%t(Zbar)
    G = fG*S
    d =  p*matrix.rank(D%*%t(D))
    TD = sqrt(d)*fG*matrix.trace(H)/matrix.trace(G) - sqrt(d)*fH

    TDsd_new = (2*(fH)*matrix.trace(S%*%S)/d-matrix.trace(S)^2/(fG*d))/(matrix.trace(S)^2/d^2)*(1+fH/fG)*fG^2/(fG^2+fG-2)
    #TDsd_old = (2*(fH)*matrix.trace(S%*%S)/d-matrix.trace(S)^2/(fG*d))/(matrix.trace(S)^2/d^2)*fG^2/(fG^2+fG-2)

    # test based on Srivastava variance estimator
    TDsrv_new = TD/sqrt(TDsd_new)
    #TDsrv_old = TD/sqrt(TDsd_old)

    # test based on U-statistics estimator
    TDsdN <- 0
    for(i in 1:a) {
      TDsdN <- TDsdN + 1/n[i]^2*Q[i, 2]
      j <- i + 1
      while(j <= a) {
        TDsdN <- TDsdN + 2*1/n[i]*1/n[j]*(matrix.trace(Si[[i]]%*%Si[[j]]))
        j <- j + 1
      }
    }
    TDsdN <- sqrt(2*fH*TDsdN*d)

    TDsdD <- 0
    for(i in 1:a) {
      TDsdD <- TDsdD + 1/n[i]^2*Q[i, 1]
      j <- i + 1
      while(j <= a) {
        TDsdD <- TDsdD + 2*1/n[i]*1/n[j]*(matrix.trace(Si[[i]])*matrix.trace(Si[[j]]))
        j <- j + 1
      }
    }
    TDsdD <- sqrt(TDsdD)
    TDsd2 <- TDsdN/TDsdD*sqrt(1+fH/fG)

    TDu = TD/TDsd2

    # Boik (1988), Rao F Approximation
    WL <- FALSE
    Lambda = abs(1/det(I(dim(H)[1])+H%*%solve(G,tol=1e-30)))

    m1 = fH
    m2 = fG
    if( (d < fG) & ((d^2+m1^2-5) > 0) & (d<=min(n)) ){

      a = sqrt((d^2+m1^2-5)/((m1*d)^2-4))
      v1 = m1*d
      v2 = a^(-1)*(m2-1/2*(d-m1+1))-1/2*(m1*d-2)
      # finite sample approximation with an F(v1, v2) distribution
      WL_test = (Lambda^(-a)-1)*v2/v1
      WL <- TRUE*(v2 > 0)
    }

    result <- data.frame(method = 1:2, fH = 1:2, fG = 1:2, test = 1:2, pvalue = 1:2)
    result$method <- c("TDsrv_new", "TDu")
    result$fH <- rep(fH, 2)
    result$fG <- rep(fG, 2)
    result$test <- c(TDsrv_new, TDu)
    result$pvalue <- as.double(sapply(result$test, FUN = function(s) { 2*min(pnorm(s), 1-pnorm(s)) }))
    result$sign.code <- sapply(as.double(result$pvalue), FUN = .hrm.sigcode )

    if(WL) {
      result[3, ] <- c("WL", v1, v2, WL_test, 1-pf(WL_test,v1,v2), .hrm.sigcode(WL_test) )
    } else {
      result[3, ] <- c("WL", NA, NA, NA, NA, NA )
    }

    return(result)
}
