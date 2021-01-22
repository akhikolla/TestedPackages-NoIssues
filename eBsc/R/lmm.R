lmm <-
function (y, q, parameters){
    neBsc <- length(y)
    basis <- BasiseBsc[[q]]
    col.ones <- rep(1,neBsc)
    DATA <- cbind(index=seq(1,length.out=neBsc),y=as.numeric(y))
    X <- basis$eigenvectorsQR[,1:q]
    Z <- basis$eigenvectorsQR[,(q+1):neBsc]%*%diag(1/sqrt(basis$eigenvalues[(q+1):neBsc]))
    N <- basis$eigenvectorsQR
    D <- diag(basis$eigenvalues)
    if(is.null(parameters)){
        lme.fit <- lme(y ~ -1+X, random = list(col.ones=(pdIdent(~-1+Z))),correlation=NULL)
        sigma2.e<- lme.fit$sigma^2
        sigma2.b<-(lme.fit$sigma^2)*exp(2*unlist(lme.fit$modelStruct))
        lambda<- ((sigma2.e)/(sigma2.b))[1] 
        V <- diag(neBsc)
        
    }else{
        lme.fit <- lme(y ~ -1+X, random = list(col.ones=(pdIdent(~-1+Z))),correlation=eval(parse(text=parameters)))
        sigma2.e <- lme.fit$sigma^2
        sigma2.b<-(lme.fit$sigma^2)*exp(2*unlist(lme.fit$modelStruct))
        lambda <- ((sigma2.e)/(sigma2.b))[1]
        corcoefs <- coef(lme.fit$modelStruct$corStruct,unconstrained=FALSE)
        V <- mycorMatrix(corcoefs,DATA)
    }
    
    list(sigma2.hat=sigma2.e,lambda.hat=log(lambda/neBsc,10),R.hat=V,X=X,Z=Z,N=N,D=D)
}
