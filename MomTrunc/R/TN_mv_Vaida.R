###############################################################################################
###############################################################################################
#This function is optimized for the case when all intervalar lower/upper censoring limits
#are finite
###############################################################################################
###############################################################################################

Vaida.IC = function(a=-c(3,2),b=c(1,2),mu=c(0,0),Sigma=matrix(c(1,-0.5,-0.5,1),2,2)){
  p=length(mu)
  if(p==1){
    s = sqrt(Sigma)
    a1 = (a-mu)/s
    b1 = (b-mu)/s
    L = pnorm(b1)-pnorm(a1)
    muY = mu+(dnorm(a1)-dnorm(b1))/L*s
    varY = Sigma+(mu-muY)*muY+(a*dnorm(a1)-b*dnorm(b1))/L*s
    return(list(mean = muY,EYY = varY+muY^2,varcov = varY))
  }else{

    #####
    Sigma = sym.matrix(Sigma)
    #####

    a1 <- (a-mu)/sqrt(diag(Sigma))
    b1 <- (b-mu)/sqrt(diag(Sigma))
    R <-  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    alpha <- pmvn.genz(lower = as.vector(a1),upper=as.vector(b1),sigma=R,uselog2 = TRUE)$Estimation
    if(2^alpha < 1e-250){
      #print("IC.Vaida corrector applied \n")
      return(corrector(a,b,mu,Sigma,bw=36))
    }
    qq = qfun.noinf(a1,b1,Sigma = R)
    da = qq$qa
    db = qq$qb
    EX <- -R%*%log2ratio(da - db,alpha)
    #EX = -R%*%(da - db)/alpha
    Eycens <- -diag(sqrt(diag(Sigma)))%*%EX+mu
    #-qq is b standardized
    if(max(abs(EX))> 10*max(abs(c(a1,b1)[is.finite(c(a1,b1))])) | any(Eycens < a | Eycens > b)){
      return(corrector(a,b,mu,Sigma,bw=36))
    }

    H =  RH = matrix(0,p,p)
    if(p==2){
      H[1,2] = H[2,1] = dmvnorm(x = as.vector(a1),sigma=R) - dmvnorm(x = as.vector(c(a1[1],b1[2])),sigma=R) - dmvnorm(x = as.vector(c(b1[1],a1[2])),sigma=R) + dmvnorm(x = as.vector(b1),sigma=R)
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }else{
      for(s in 1:(p-1)){
        for(t in (s+1):p){
          ##print(s);#print(t)
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          Sj = R[-c(s,t), c(s,t), drop=F]%*%invR
          V  =  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] = H[t,s] = log2prod0(dmvnorm(x = a1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                                      pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%a1[c(s,t),drop=F]),sigma=V,uselog2 = TRUE)$Estimation) -
            log2prod0(dmvnorm(x = c(a1[s],b1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                      pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(a1[s],b1[t])),sigma=V,uselog2 = TRUE)$Estimation) -
            log2prod0(dmvnorm(x = c(b1[s],a1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                      pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(b1[s],a1[t])),sigma=V,uselog2 = TRUE)$Estimation) +
            log2prod0(dmvnorm(x = b1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                      pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%b1[c(s,t),drop=F]),sigma=V,uselog2 = TRUE)$Estimation)
          RH[s,t] = RH[t,s] = R[s,t]*H[s,t]
        }
      }
    }
    a1[is.infinite(a1)] = 0
    b1[is.infinite(b1)] = 0
    h = (a1*da - b1*db) - apply(RH, 1, sum)
    diag(H) = h
    #EX <- -R%*%(da - db)/alpha   # a vector with a length of p
    EXX <- R + R%*%log2ratio(H,alpha)%*%R
    varX <- EXX-EX%*%t(EX)
    varyic <- diag(sqrt(diag(Sigma)))%*%varX%*%diag(sqrt(diag(Sigma)))
    E2yy <- varyic+Eycens%*%t(Eycens)

    #Validating positive variances
    bool = diag(varyic) < 0
    if(sum(bool)>0){
      #print("negative variance found")
      out = corrector(a,b,mu,Sigma,bw=36)
      out$mean = Eycens
      out$EYY = out$varcov + out$mean%*%t(out$mean)
      return(out)
    }

  }
  return(list(mean=Eycens,EYY=E2yy,varcov=varyic))
}

###############################################################################################
###############################################################################################
#This function is optimized for the case when it DOES exist infinite values in the lower/upper
#truncation limits
###############################################################################################
###############################################################################################

Vaida.LRIC<-function(a=-c(-Inf,2),b=c(1,Inf),mu=c(0,0),Sigma=matrix(c(1,-0.5,-0.5,1),2,2)){
  p=length(mu)
  if(p==1){
    s = sqrt(Sigma)
    a1 = (a-mu)/s
    b1 = (b-mu)/s
    L = pnorm(b1)-pnorm(a1)
    muY = mu+(dnorm(a1)-dnorm(b1))/L*s
    if(a == -Inf){
      a = 0
    }
    if(b == Inf){
      b = 0
    }
    varY = Sigma+(mu-muY)*muY+(a*dnorm(a1)-b*dnorm(b1))/L*s
    return(list(mean = muY,EYY = varY+muY^2,varcov = varY))
  }else{

    #####
    Sigma = sym.matrix(Sigma)
    #####

    a1 <- (a-mu)/sqrt(diag(Sigma))
    b1 <- (b-mu)/sqrt(diag(Sigma))
    R <-  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))


    alpha <- pmvn.genz(lower = as.vector(a1),upper=as.vector(b1),sigma=R,uselog2 = TRUE)$Estimation
    if(2^alpha < 1e-250){
      #print("IC.Vaida corrector applied \n")
      return(corrector(a,b,mu,Sigma,bw=36))
    }
    qq = qfun(a1,b1,Sigma = R)
    da = qq$qa
    db = qq$qb
    EX <- -R%*%log2ratio(da - db,alpha)
    Eycens <- -diag(sqrt(diag(Sigma)))%*%EX+mu

    #-qq is b standardized
    if(max(abs(EX))> 10*max(abs(c(a1,b1)[is.finite(c(a1,b1))])) | any(Eycens < a | Eycens > b)){
      return(corrector(lower = a,upper = b,mu,Sigma,bw=36))
    }

    #var
    H =  RH = matrix(0,p,p)
    if(p==2){
      H[1,2] = H[2,1] = dmvnorm(x = as.vector(a1),sigma=R) - dmvnorm(x = as.vector(c(a1[1],b1[2])),sigma=R) - dmvnorm(x = as.vector(c(b1[1],a1[2])),sigma=R) + dmvnorm(x = as.vector(b1),sigma=R)
      #sigma==R since b1 is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }else{
      for(s in 1:(p-1)){
        for(t in (s+1):p){
          ##print(s);#print(t)
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          Sj = R[-c(s,t), c(s,t), drop=F]%*%invR
          V  =  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] = H[t,s] = ifelse(any(is.infinite(a1[c(s,t)])),0,
                                   log2prod0(dmvnorm(x = a1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                                             pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%a1[c(s,t),drop=F]),sigma=V,uselog2 = TRUE)$Estimation)
          ) -
            ifelse(any(is.infinite(c(a1[s],b1[t]))),0,
                   log2prod0(dmvnorm(x = c(a1[s],b1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                             pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(a1[s],b1[t])),sigma=V,uselog2 = TRUE)$Estimation)
            ) -
            ifelse(any(is.infinite(c(b1[s],a1[t]))),0,
                   log2prod0(dmvnorm(x = c(b1[s],a1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                             pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(b1[s],a1[t])),sigma=V,uselog2 = TRUE)$Estimation)
            ) +
            ifelse(any(is.infinite(b1[c(s,t)])),0,
                   log2prod0(dmvnorm(x = b1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),
                             pmvn.genz(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%b1[c(s,t),drop=F]),sigma=V,uselog2 = TRUE)$Estimation)
            )
          RH[s,t] = RH[t,s] = R[s,t]*H[s,t]
        }
      }
    }
    a1[is.infinite(a1)] = 0
    b1[is.infinite(b1)] = 0
    h = (a1*da - b1*db) - apply(RH, 1, sum)
    diag(H) = h
    #EX <- -R%*%(da - db)/alpha   # a vector with a length of p
    EXX <- R + R%*%log2ratio(H,alpha)%*%R
    #EXX <- R + R%*%(H/alpha)%*%R
    varX <- EXX-EX%*%t(EX)
    #Eycens <- -diag(sqrt(diag(Sigma)))%*%EX+mu
    varyic <- diag(sqrt(diag(Sigma)))%*%varX%*%diag(sqrt(diag(Sigma)))
    E2yy <- varyic+Eycens%*%t(Eycens)

    #Validating positive variances
    bool = diag(varyic) < 0
    if(sum(bool)>0){
      #print("negative variance found")
      out = corrector(a,b,mu,Sigma,bw=36)
      out$mean = Eycens
      out$EYY = out$varcov + out$mean%*%t(out$mean)
      return(out)
    }
  }
  return(list(mean=Eycens,EYY=E2yy,varcov=varyic))
}


###############################################################################################
###############################################################################################
#This is the original Vaida's function for upper truncation (right censoring)
###############################################################################################
###############################################################################################

Vaida.RC = function(b=c(1,2),mu=c(0,0), Sigma = diag(2)){
  p=length(mu)
  if (p==1) {
    qq = (1/sqrt(Sigma))*(-b+mu)
    R = 1
    alpha = pnorm(-qq)
    dd = dnorm(-qq)
    H = qq*dd
    EX = dd/alpha   # a vector with a length of p
    EXX = 1+H/alpha
    varX = EXX-EX^2
    Eycens = -sqrt(Sigma)*EX+mu
    varyic= varX*Sigma
    E2yy=varyic+Eycens^2
  }
  else {
    #####
    Sigma = sym.matrix(Sigma)
    #####

    qq = diag(1/sqrt(diag(Sigma)))%*%(-b+mu)
    R =  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    alpha = pmvn.genz(upper=as.vector(-qq),sigma=R,uselog2 = TRUE)$Estimation
    if(2^alpha < 1e-250){
      #print("RC.Vaida corrector applied \n")
      return(corrector(upper = b,mu = mu,Sigma = Sigma,bw=36))
    }
    ##print(qq)
    dd = rep(0, p)   #derivative vector
    for (j in 1:p){
      V = R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu = -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] = log2prod0(dnorm(-qq[j]),tlrmvnmvt::pmvn(upper=as.vector(nu),sigma=V,uselog2 = TRUE))
    }
    EX = R%*%log2ratio(dd,alpha)
    Eycens = -diag(sqrt(diag(Sigma)))%*%EX+mu
    #-qq is b standardized
    if(max(abs(EX))> 10*max(abs(qq[is.finite(qq)])) | any(Eycens > b)){
      return(corrector(upper = b,mu = mu,Sigma = Sigma,bw=36))
    }

    H = matrix(rep(0, p*p), nrow=p)
    RH = matrix(rep(0, p*p), nrow=p)
    if(p==2){
      H[1,2] = H[2,1] = dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] = RH[2,1] = R[1,2]*H[1,2]
    }else{
      for( s in 1:(p-1)){
        for (t in (s+1):p){
          invR = solve(R[c(s,t), c(s,t), drop=F])
          nu = -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V =  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] = H[t,s] = log2prod0(dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2)),tlrmvnmvt::pmvn(upper=as.vector(nu), sigma=V,uselog2 = TRUE))
          RH[s,t] = RH[t,s] = R[s,t]*H[s,t]
        }
      }
    }
    h = qq*dd-apply(RH, 1, sum)
    diag(H) = h
    #EX = R%*%dd/alpha
    EXX <- R + R%*%log2ratio(H,alpha)%*%R
    #EXX = R + R%*%(H/alpha)%*%R
    varX = EXX-EX%*%t(EX)
    varyic = diag(sqrt(diag(Sigma)))%*%varX%*%diag(sqrt(diag(Sigma)))
    E2yy = varyic+Eycens%*%t(Eycens)

    #Validating positive variances
    bool = diag(varyic) < 0
    if(sum(bool)>0){
      #print("negative variance found")
      out = corrector(upper = b,mu = mu,Sigma = Sigma,bw=36)
      out$mean = Eycens
      out$EYY = out$varcov + out$mean%*%t(out$mean)
      return(out)
    }
  }
  return(list(mean=Eycens,EYY=E2yy,varcov=varyic))
}
