EBCparallel<-function(y, qq = NULL, method = "N", parallel = FALSE,  R0 = NULL, zero_range = c(-45,1), ARpMAq = c(0,0), trace = FALSE, tol.lambda = 0.01, tol.rho = 0.01, max.iter = 50){                    

    ##constants
    neBsc <-length(y)
    EeBsc <- exp(1)
    
    ##initialize variables
    f.Hat<- R.Hat <- lambda.Hat <- sigma2.Hat <- niter <- NA

    ##one q

    if(!is.null(qq)){
        ##basis
        useful <- FALSE
        if(exists("BasiseBsc")){if(length(BasiseBsc)==max(qq,2)){ch <- sample(min(qq,2),1); useful = identical(BasiseBsc[[ch]],drbasis(nn=neBsc,qq=ch))}}
        if(!useful){
            BasiseBsc<-list()
            for (i in union(2,qq)){BasiseBsc[[i]]=drbasis(nn=neBsc,qq=i)}
            BasiseBsc <<- BasiseBsc
        }

    
        ##computation
        O <- EBCq(y = y, q = qq, method = method, R0 = R0, zero_range = zero_range, ARpMAq = ARpMAq, tol.lambda = tol.lambda, tol.rho = tol.rho, max.iter = max.iter)
        f.Hat <- O$f.hat
        R.Hat <- O$R.hat
        lambda.Hat <- O$lambda.hat
        sigma2.Hat <- O$sigma2.hat        
        etq.Hat <- O$etq.hat
    
        q.Hat <- qq
        niter <- O$niter

        ##message
        if(trace & niter<50){
            message(paste("eBsc success!: convergence achieved after",niter,"iterations."),"\n",
                paste("eBsc warning!: estimation done for the provided q = ",qq,".",sep=""),"\n")}
        if(trace & niter==50){
            message(paste("eBsc fail!: no convergence achieved after",niter,"iterations."),"\n",
                paste("eBsc fail!: saved results obtained at iteration 50."),"\n")}    
    }

    ##all q

    if(is.null(qq)){

        ##basis
        useful <- FALSE
        if(exists("BasiseBsc")){if(length(BasiseBsc)==6){ch <- sample(6,1); useful = identical(BasiseBsc[[ch]],drbasis(nn=neBsc,qq=ch))}}
        if(!useful){
            BasiseBsc <<- drbasis(nn=neBsc)
            if(trace){message("eBsc success!: DR basis computed.","\n")}        
        }
        
        ##computation
        if(parallel){                        
            nworkers <- detectCores()
            cl <- makePSOCKcluster(nworkers, varlist = ls())
            clusterExport(cl, c("neBsc","EeBsc","BasiseBsc","get_lambda","smoothrhos","Eq"))
            O <- parApply(cl = cl, matrix(1:6, nrow=1), 2, EBCq, y = y, method = method, R0 = R0, zero_range = zero_range, ARpMAq = ARpMAq,
                                  tol.lambda = tol.lambda, tol.rho = tol.rho, max.iter = max.iter)
            stopCluster(cl)
        }else{
            O <- list()
            for(j in 1:6){O[[j]] <- EBCq(y = y, q = j, method = method, R0 = R0, zero_range = zero_range, ARpMAq = ARpMAq,
                                         tol.lambda = tol.lambda, tol.rho = tol.rho, max.iter = max.iter)}
        }
    
        etq.Hat <- NULL; for(j in 1:6){etq.Hat <- c(etq.Hat,O[[j]]$etq.hat)}
        q.Hat <- get_q(etq.Hat);    

        if(q.Hat!=0){
            f.Hat      <- O[[q.Hat]]$f.hat
            R.Hat      <- O[[q.Hat]]$R.hat
            lambda.Hat <- O[[q.Hat]]$lambda.hat
            sigma2.Hat <- O[[q.Hat]]$sigma2.hat
            niter      <- O[[q.Hat]]$niter
        }else{
            message(" eBsc fail!: cannot select optimal q.","\n",
                paste("eBsc fail!: reporting results for each q in << out >> object."),"\n");
        }

        ##message
        if(trace & niter<50){
            message(paste("eBsc success!: convergence achieved after",niter,"iterations."),"\n")}
        if(trace & niter==50){
            message(paste(" eBsc fail!: no convergence achieved after",niter,"iterations."),"\n",
                paste("eBsc fail!: saved results obtained at iteration 50"),"\n")}    
    }

    return(list(
        f.hat = f.Hat,
        R.hat = R.Hat,    
        lambda.hat = lambda.Hat,
        sigma2.hat =sigma2.Hat,        
        etq.hat = etq.Hat,
        q.hat = q.Hat,
        niter = niter,
        data = y,
        res = y-f.Hat,
        out = O))
}
