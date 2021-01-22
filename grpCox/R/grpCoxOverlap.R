#############################################
#####  Overlapping Group Cox  function  #####
#############################################
grpCoxOverlap <- function(X0, y, group, penalty=c("glasso", "gSCAD", "gMCP"),lambda=NULL, nlambda=100, rlambda=NULL,
                         gamma=switch(penalty, gSCAD = 3.7, 3),standardize = TRUE, thresh=1e-3, maxit=1e+4,
                         returnLatent = TRUE){
  X <- expandedlatent(X0, group)
  
  g <- attr(X,'grp_vec')
  m <- sqrt(table(g[g!=0]))
  dt <- sglcoxinfo(X,y,g,m,standardize = standardize)
  
  N0 <- dim(dt$x)[1]
  P <- dim(dt$x)[2]
  
  penalty <- match.arg(penalty)
  # checking
  if (nlambda < 2) stop("nlambda must be at least 2")
  
  #### lambda path
  if(is.null(lambda)){
    lambda_max  <-  max_lambda(dt)
    if(is.null(rlambda)){
      rlambda <- ifelse(N0>P, 0.001, 0.05)
    }
    lambda_min <- lambda_max*rlambda
    lambda <- lambda_max*(lambda_min/lambda_max)^(c(0:(nlambda-1))/(nlambda-1))
  }else{
    nlambda <- length(lambda)
  }
  
  #### Main function
  out  <-  grpCoxQ(dt, penalty, lambda, nlambda, gamma, thresh, maxit)
  
  out_betas <- matrix(rep(0,nlambda*P), nrow=P)
  ind <- intersect(intersect(which(out$flag==0), which(out$ll != "NaN")), which(out$ll != 0))
  lambdai <- lambda[ind]
  nlambdai <- length(ind)
  if(nlambdai==0) return(NULL)
  out$betaSTD <- out$beta[,ind]
  ## Check the group order 
  if(is.null(dt$ord.inv)){
    ## Get the beta without standardizing
    for(i in ind){
      beta0_i <- out$beta[,i]
      xscale_q <- dt$scale
      for(p in (1:P)){
        beta0_i[p] <- beta0_i[p]/xscale_q[p] 
      }
      out_betas[,i] <- beta0_i
    }
    out$beta_O <- out_betas[,ind]
  }else{
    ## Get the beta without standardizing
    for(i in ind){
      beta0_i <- out$beta[,i]
      xscale_q <- dt$scale
      for(p in (1:P)){
        beta0_i[p] <- beta0_i[p]/xscale_q[p] 
      }
      out_betas[,i] <- beta0_i[dt$ord.inv]
      out$beta[,i]  <-  out$beta[,i][dt$ord.inv]
    }
    out$beta_O <- out_betas[,ind]
  }
  out$betaOri <- latent2ori(coef.latent=out$beta_O, var_grp=attr(X,'var_grp'),grp_vec=attr(X,'grp_vec'))
  
  if(returnLatent){
    return(list(aBetaLatent=out$beta_O, aBetaOri=out$betaOri, lambda=lambdai, 
                ll=out$ll[ind], group=group, glatent=g))
  }else{
    return(list(aBetaOri=out$betaOri, lambda=lambdai, 
                ll=out$ll[ind], group=group, glatent=g))
  }
}

#####################################################
#####  Cross-validation Overlapping Group Cox  ######
#####################################################
cv.grpCoxOverlap <- function(X0, y, group, penalty=c("glasso", "gSCAD", "gMCP"),lambda=NULL, nlambda=100,
                            rlambda=NULL,gamma=switch(penalty, SCAD = 3.7, 3), 
                            standardize = TRUE, thresh=1e-3, maxit=1e+4,
                            nfolds=10, foldid=NULL, returnLatent = TRUE){
  X <- expandedlatent(X0, group)
  
  g <- attr(X,'grp_vec')
  m <- sqrt(table(g[g!=0]))
  dt <- sglcoxinfo(X,y,g,m,standardize = standardize)
  
  N0 <- dim(dt$x)[1]
  P <- dim(dt$x)[2]
  
  penalty <- match.arg(penalty)
  # checking
  if (nlambda < 2) stop("nlambda must be at least 2")
  
  #### lambda path
  if(is.null(lambda)){
    lambda_max  <-  max_lambda(dt)
    if(is.null(rlambda)){
      rlambda <- ifelse(N0>P, 0.001, 0.05)
    }
    lambda_min <- lambda_max*rlambda
    lambda <- lambda_max*(lambda_min/lambda_max)^(c(0:(nlambda-1))/(nlambda-1))
  }else{
    nlambda <- length(lambda)
  }
  
  #### Main function
  out  <-  grpCoxQ(dt, penalty, lambda, nlambda, gamma, thresh, maxit)
  
  out_betas <- matrix(rep(0,nlambda*P), nrow=P)
  ind <- intersect(intersect(which(out$flag==0), which(out$ll != "NaN")), which(out$ll != 0))
  lambdai <- lambda[ind]
  nlambdai <- length(ind)
  if(nlambdai==0) return(NULL)
  out$betaSTD <- out$beta[,ind]
  ## Check the group order 
  if(is.null(dt$ord.inv)){
    ## Get the beta without standardizing
    for(i in ind){
      beta0_i <- out$beta[,i]
      xscale_q <- dt$scale
      for(p in (1:P)){
        beta0_i[p] <- beta0_i[p]/xscale_q[p] 
      }
      out_betas[,i] <- beta0_i
    }
    out$beta_O <- out_betas[,ind]
  }else{
    ## Get the beta without standardizing
    for(i in ind){
      beta0_i <- out$beta[,i]
      xscale_q <- dt$scale
      for(p in (1:P)){
        beta0_i[p] <- beta0_i[p]/xscale_q[p] 
      }
      out_betas[,i] <- beta0_i[dt$ord.inv]
      out$beta[,i]  <-  out$beta[,i][dt$ord.inv]
    }
    out$beta_O <- out_betas[,ind]
  }
  out$betaOri <- latent2ori(coef.latent=out$beta_O, var_grp=attr(X,'var_grp'),grp_vec=attr(X,'grp_vec'))
  
  ## Check whether or not to do cross-validation
  if(nfolds==1 & is.null(foldid)){
    fit <- data.frame(lambda=lambdai)
    return(list(aBetaSTD=out$betaSTD, aBetaO=out$beta_O, fit=fit, lambda=lambdai, iter=out$iter, nz=out$nz))
  }else{
    if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")
    #### Split data for cross-validation
    foldid<-sample(rep(seq(nfolds), length=N0))
    n<-seq(1:N0)
    ### Do cross-validation
    cvPL <- matrix(NA, nrow=nfolds, ncol=nlambdai)
    outi <- list()
    for(i in 1:nfolds){
      ni <- n[which(foldid!=i)]
      Xi <- X[ni,]; yi = y[ni,]
      di <- sglcoxinfo(Xi,yi,g,m,standardize = standardize)
      outi[[i]] <- cvgrpCoxQ(di, penalty, lambdai, nlambdai, gamma, thresh, maxit, dt)
      cvPL[i, 1:outi[[i]]$nlambda] <- outi[[i]]$lf[1:outi[[i]]$nlambda] - outi[[i]]$ll[1:outi[[i]]$nlambda]
    }
    
    #### Process results
    cvPL <- matrix(cvPL[,1:nlambdai], ncol = nlambdai)
    cvraw <- cvPL; nfoldi <- apply(!is.na(cvraw), 2, sum); rm(cvPL)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvse <- sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, mean, na.rm = TRUE)/(nfoldi-1))
    
    indexi <- which.max(cvm)
    
    ######penalized extensions
    l <- rev(lambdai)
    cvl <- rev(cvm*N0)
    #lambda maximizes cvl and its corresponding cvl, number of non-zero coefs
    l.cvl <- lambdai[indexi]
    cvlMax <- cvl[which(l==l.cvl)]
    #degree of freedom = number of non-zeros
    df <- rep(0,nlambdai)
    for(j in (1:nlambdai)){
      df[j] <- sum(out$betaSTD[,j] != 0)  
    }
    nz <- rev(df)
    nzMax <- nz[which(l==l.cvl)]
    
    nzmin <- min(nz)
    cvl0 <- max(cvl[which(nz==nzmin)])
    
    if(nzMax != 0){
      #lasso-pcvl: penalized cross-validated log-likelihood
      pcvl <- cvl - (((cvlMax - cvl0)/nzMax) *nz)
      l.pcvl <- l[which.max(pcvl)]
    }else{
      l.pcvl <- l.cvl
    }
    
    indexm <- which(lambdai==l.pcvl)
    CV.max <- cvm[indexm]
    
    
    #################### Result
    temi <- rep("", nlambdai)
    temi[indexi] <- "cvmax"
    if(indexm==indexi){
      temi[indexm] <- "pcvl=cvmax"
    }else{
      temi[indexm] <- "pcvl"
    }
    
    temCV <- data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, index=temi, stringsAsFactors=FALSE)
    
    if(returnLatent){
      return(list(aBetaLatent=out$beta_O, aBetaOri=out$betaOri, pBetaLatent=out$beta_O[,indexm],pBetaOri=out$betaOri[,indexm], 
                  mBetaLatent=out$beta_O[,indexi],mBetaOri=out$betaOri[,indexi],fit=temCV, lambda=lambdai,  
                  group=group,glatent=g,cvmax=CV.max,lambda.max=lambdai[indexi], lambda.pcvl=lambdai[indexm]))
    }else{
      return(list(aBetaOri=out$betaOri, pBetaOri=out$betaOri[,indexm], 
                  mBetaOri=out$betaOri[,indexi],fit=temCV, lambda=lambdai, 
                  group=group,glatent=g,cvmax=CV.max,lambda.max=lambdai[indexi], lambda.pcvl=lambdai[indexm]))
    }
  }
}