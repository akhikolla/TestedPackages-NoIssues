EBCq <- function(y, q, method = "N", R0 = NULL, zero_range = c(-45,1), ARpMAq = c(0,0), tol.lambda = 0.01, tol.rho = 0.01, max.iter = 50)
{
    ##input
    neBsc <- length(y)
    sdy <- sd(y)
    y <- y/sdy
    
    y.vec <- matrix(y, ncol = 1)
    x <- seq(1, neBsc, length = neBsc)
    Phi.mat <- BasiseBsc[[q]]$eigenvectorsQR
    eta.vec <- as.vector(BasiseBsc[[q]]$eigenvalues)
    coefs <- t(Phi.mat) %*% y.vec
    eigens <- neBsc*eta.vec 
    
    ##default
    Rseq <- array(data = NA, dim = c(25+1,neBsc,neBsc));
    if(is.null(R0)){R0 <- diag(neBsc)}
    Rseq[1,,] <- R0;
    

    ##constants
    eps <- .Machine$double.eps

    ##initialize variables
    lambda.hat <- sigma2.hat <- etq.hat <- NA
    f.hat <- R.eigen.values <- rep(NA,neBsc);  
    R.eigen.values.seq <- matrix(NA,neBsc,25)
    log10lambdas <- vars <- lambdas.hat <- sigma2s.hat <- rep(NA,25)
    S.hat <- R.hat <- matrix(NA,neBsc,neBsc)
    niter <- NA    
   
#####################
#deterministic method
#####################
    
    if(method == "D"){
        
        R.hat <- R0
        
        if (all.equal(R.hat,diag(neBsc))){
            rho.fit = 1
        }else{
            s = 1:(neBsc-1);
            rho.fit = 1+2*sapply(1:neBsc, function(i) sum(cos(pi*(i-1)*s/(neBsc-1))*R.hat[1,2:neBsc]))
        }

        log10lambda.hat <- get_lambda(zero_range = zero_range, coefs = coefs, eigens = eigens, rhos = rho.fit , neBsc = neBsc, q = q)        
        Lambda = 10^log10lambda.hat

        ##raw eigenvalues
        rho = coefs ^ 2 * Lambda * eigens * rho.fit / (1 + Lambda * eigens * rho.fit)
        #rho[rho > quantile(rho, 0.995)] = quantile(rho, 0.995)
        rho = rho / mean(rho)
      
        ##smoothed rhos
        rho.fit1 = smoothrhos(rho = rho, BasiseBsc = BasiseBsc[[2]], zero_range = zero_range)        
        rho.fit = rho.fit1/mean(rho.fit1)
        
        ##final estimates
        f.hat   = Phi.mat %*% (coefs/(1+Lambda * eigens * rho.fit))
        etq.hat = Eq(log10lambda = log10(Lambda), coefs = coefs, eigens = eigens, rhos = rho.fit, neBsc = neBsc, q = q)
        sigma2.hat  = ( sum( coefs ^ 2 * Lambda * eigens / (1 + Lambda * eigens * rho.fit)) + 1)/( neBsc + 1)
          
        outcome <-
            list(
                f.hat = f.hat*sdy,
                R.hat = R.hat[1,],
                lambda.hat = Lambda,
                etq.hat = etq.hat,
                niter = 1,
                sigma2.hat = (sdy^2)*sigma2.hat
            )
    }
    
#####################
#nonparametric method
#####################

    if(method=="N"){
        R.hat <- R0
        
        if (all.equal(R.hat,diag(neBsc))){
            rho.fit = 1
        }else{
            s = 1:(neBsc-1);
            rho.fit = 1+2*sapply(1:neBsc, function(i) sum(cos(pi*(i-1)*s/(neBsc-1))*R.hat[1,2:neBsc]))
        }

        log10lambda.hat <- get_lambda(zero_range = zero_range, coefs = coefs, eigens = eigens, rhos = rho.fit , neBsc = neBsc, q = q)        
        Lambda = 10^log10lambda.hat
        
        count = 1
        
        repeat
        {
            ##Raw eigenvalues
            rho = coefs ^ 2 * Lambda * eigens * rho.fit / (1 + Lambda * eigens * rho.fit)
            #rho[rho > quantile(rho, 0.995)] = quantile(rho, 0.995)
            rho = rho / mean(rho)
      
            ##Smoothed rhos
            rho.fit1 = smoothrhos(rho = rho, BasiseBsc = BasiseBsc[[2]], zero_range = zero_range)
            rho.fit1 = rho.fit1/mean(rho.fit1)
            
            ##New Lambda
            log10lambda.hat = get_lambda(zero_range = zero_range, coefs = coefs, eigens = eigens, rhos = rho.fit1 , neBsc = neBsc, q = q)        
            Lambda1         = 10^log10lambda.hat
      
            ##check convergence of both estimating equaiton (for rho and lambda)
            Lambda.diff = abs(log(Lambda1,neBsc) - log(Lambda,neBsc))
            rho.diff    = mean(abs(rho.fit - rho.fit1))
      
            if (((Lambda.diff<tol.lambda) & (rho.diff<tol.rho))|count==max.iter ){break}
      
            count = count + 1
        }

        Lambda  = Lambda1
        rho.fit = rho.fit1

        
        ##final estimates
        f.hat   = Phi.mat %*% (coefs / (1 + Lambda * eigens * rho.fit))
        etq.hat = Eq(log10lambda = log10(Lambda), coefs = coefs, eigens = eigens, rhos = rho.fit, neBsc = neBsc, q = q)
        R.hat   = toeplitz(sapply(1:neBsc, function(i) mean(cos(pi * seq(0,1,length.out=neBsc) * (i-1)) * rho.fit)))
        sigma2.hat  = (sum(coefs ^ 2 * Lambda * eigens/(1 + Lambda * eigens * rho.fit)) + 1)/(neBsc + 1)
        
        outcome<-
            list(
                f.hat = sdy*f.hat,
                R.hat = R.hat[1,],
                lambda.hat = Lambda,
                etq.hat = etq.hat,
                niter = count,
                sigma2.hat = (sdy^2)*sigma2.hat
            )
    }


##################
#parametric method
##################

if(method=="P"){
    if((ARpMAq[1] == 0)&(ARpMAq[2] == 0)){
        correlation <- NULL;
        mm <- lmm(y = y, q = q, parameters = correlation);
    }else{
        correlation <- parse(text = paste("corARMA(","p=",ARpMAq[1],",q=",ARpMAq[2],")",sep=""))
        mm <- lmm(y = y, q = q, parameters = correlation)
    }
    R.hat      <- mm$R.hat
    sigma2.hat <- mm$sigma2.hat

    if (norm((R.hat-diag(neBsc)),"F")<2e-06){
        rho.fit = 1
    }else{
        s = 1:(neBsc-1);
        rho.fit = 1+2*sapply(1:neBsc, function(i) sum(cos(pi*(i-1)*s/(neBsc-1))*R.hat[1,2:neBsc]))
    }

    log10lambda.hat = mm$lambda.hat
    Lambda = 10^log10lambda.hat

    rho = coefs ^ 2 * Lambda * eigens * rho.fit / (1 + Lambda * eigens * rho.fit)
    #rho[rho>quantile(rho,0.995)] = quantile(rho, 0.995) 
    rho = rho / mean(rho)
      
    ##smoothed rhos
    rho.fit1 = smoothrhos(rho = rho, BasiseBsc = BasiseBsc[[2]], zero_range = zero_range)
    rho.fit = rho.fit1/mean(rho.fit1)  
    
    ##final estimates
    f.hat   = Phi.mat %*% (coefs/(1+Lambda * eigens * rho.fit))
    etq.hat = Eq(log10lambda = log10(Lambda), coefs = coefs, eigens = eigens, rhos = rho.fit, neBsc = neBsc, q = q)
    R.hat   = sapply(1:neBsc, function(i) mean(cos(pi * x * (i-1)) * rho.fit))
            
    outcome <-
        list(
            f.hat = sdy*f.hat,
            R.hat = mm$R.hat[1,],
            lambda.hat = Lambda,
            etq.hat = etq.hat,
            niter = 1,
            sigma2.hat = (sdy^2)*sigma2.hat
        )    
}
    
outcome
}
