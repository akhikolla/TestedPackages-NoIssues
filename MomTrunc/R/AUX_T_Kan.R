Gauss2F1b <- function(a,b,z){
  if(z>=0 & z<1){
    genhypergeo(a,b,z)
  }else{
    genhypergeo(c(a[1],b-a[2]),b,1-1/(1-z))/(1-z)^a[1]
  }
}

# %
# %  This program computes the function
# %  L = Gamma(nu/2)*L_nu(a,b;0,S/nu), allowing nu
# %  to be less than or equal to zero.
# %
hfun = function(a,b,mu,S,nu){
  # print("hfun")
  # print(nu)
  n = length(a)
  s = sqrt(diag(as.matrix(S)))
  as = (a-mu)/s
  bs = (b-mu)/s
  R  = S/(s%*%t(s))
  c0 = lfun(as,bs,R,nu)
  if(n==1){
    if ((is.infinite(a)||is.infinite(b)) && nu<=1){
      y = NaN
    }else{
      y = mu*c0
      if(nu==1){
        c = s/(2*sqrt(pi))
      }else{
        c = s/(2*sqrt(pi))*gamma((nu-1)/2)
      }

      if(is.finite(a)){
        if(nu==1){
          y = y-c*log(1+as^2)
        }else{
          y = y+c*(1+as^2)^(-(nu-1)/2)
        }
      }

      if(is.finite(b)){
        if(nu==1){
          y = y+c*log(1+bs^2)
        }else{
          y = y-c*(1+bs^2)^(-(nu-1)/2)
        }
      }
    }
    return(y)
  }

  if(nu==0){
    # if(n<=4){
    #   y = matrix(NaN,n,1)
    #   for(i in 1:n){
    #     mu2 = mu[-i]
    #     a2 = a[-i]-mu2
    #     b2 = b[-i]-mu2
    #     if(is.finite(a[i]) && is.finite(b[i])){
    #       v1 = S[i,i]
    #       c = S[-i,i]/v1;
    #       R2 = S[-i,-i]-c%*%t(S[i,-i])
    #       s2 = sqrt(diag(as.matrix(R2)))
    #       R2 = R2/(s2%*%t(s2))
    #
    #       #integral
    #       integd = function(x1){
    #         #n1 = length(x1)
    #         x1a = x1-mu[i]
    #         sx1 = sqrt(1+x1a^2/v1)
    #         a2m = (a2-c*x1a)/(s2*sx1)
    #         b2m = (b2-c*x1a)/(s2*sx1)
    #         return(pmvt(lower = c(a2m),upper = c(b2m),df = 1,corr = R2)[1]*x1/sx1)
    #       }
    #
    #       integd2 = function(x2){
    #         return(sapply(x2,integd))
    #       }
    #
    #       #seqq = seq(a[i],b[i],length.out = 2000)
    #       #y[i] = (b[i] - a[i])*mean(integd2(seqq))/sqrt(v1)
    #
    #       y[i] = tryCatch(integrate(integd2,a[i],b[i],abs.tol = 1e-4,rel.tol = 1e-4)$value/sqrt(v1),
    #                       error=function(e){
    #                         seqq = seq(a[i],b[i],length.out = 2000);
    #                         y[i] = (b[i] - a[i])*mean(integd2(seqq))/sqrt(v1)})
    #
    #       # % When a2 and b2 has two finite pairs or above, we can also do integration but
    #       # % often numerical integration is not accurate, so rely on linear interpolation
    #       # % instead.
    #     }else{
    #       if(sum(is.finite(a) & is.finite(b))<2){
    #         y[i] = NaN
    #       }
    #     }
    #   }
    # }else{
    y = (hfun(a,b,mu,S,-1e-1)+hfun(a,b,mu,S,1e-1))/2
    #y = (hfun(a,b,mu,S,-1)+hfun(a,b,mu,S,1))/2
    #}
  }else{
    if(nu==1){
      qq = qfunT(a-mu,b-mu,S,1)
      y = c0*mu + sqrt(pi)*S%*%(qq$qa-qq$qb)
    }else{
      La = Lb = matrix(0,n,1)
      for(i in 1:n){
        SS = S[-i,-i]-S[-i,i]%*%t(S[i,-i])/S[i,i]
        if(!is.infinite(as[i])){
          ma = mu[-i]+S[-i,i]*as[i]/s[i]
          La[i] = (1+as[i]^2)^(-(nu-1)/2)*lfun(a[-i]-ma,b[-i]-ma,SS*(1+as[i]^2),nu-1)
        }
        if(!is.infinite(bs[i])){
          mb = mu[-i]+S[-i,i]*bs[i]/s[i]
          Lb[i] = (1+bs[i]^2)^(-(nu-1)/2)*lfun(a[-i]-mb,b[-i]-mb,SS*(1+bs[i]^2),nu-1)
        }
      }
      y = c0*mu + 1/(2*sqrt(pi))*S%*%((La-Lb)/s)
    }
  }
  return(y)
}


############################################################################################################################################
############################################################################################################################################
############################################################################################################################################


# %
# %  This program computes the function
# %  L = Gamma(nu/2)*L_nu(a,b;0,S/nu), allowing nu
# %  to be less than or equal to zero.
# %
lfun = function(a,b,S,nu){
  # print("lfun")
  # print(nu)
  n = length(a)
  s = sqrt(diag(as.matrix(S)))
  as = a/s
  bs = b/s
  if(n==1){
    if(nu>0){
      y = gamma(nu/2)*(pt(sqrt(nu)*bs,nu)-pt(sqrt(nu)*as,nu))
    }else{
      if(nu==0){
        y = asinh(bs)-asinh(as)
      }else{
        if (nu==-1){
          y = NaN
        }else{
          if(nu>-2){
            y = gamma(nu/2)*(pt(sqrt(nu+2)*bs,nu+2)-pt(sqrt(nu+2)*as,nu+2))
            cc = gamma((nu+1)/2)/(nu*sqrt(pi))
            if(is.finite(a)){y = y+cc*as/(1+as^2)^((nu+1)/2)}
            if(is.finite(b)){y = y-cc*bs/(1+bs^2)^((nu+1)/2)}
          }else{
            y = gamma((nu+1)/2)/sqrt(pi)*(bs*Gauss2F1b(c(1/2,(nu+1)/2),3/2,-bs^2)-as*Gauss2F1b(c(1/2,(nu+1)/2),3/2,-as^2))
          }
        }
      }
    }
    return(y)
  }
  #  When some of the pairs of (a(i),b(i))=(-Inf,Inf), we can
  #  reduce the problem to a lower dimension.
  ind = !(is.infinite(as) & is.infinite(bs))
  #at least one is finite
  if(sum(ind)<n){
    if(sum(ind)==0){
      y = gamma(nu/2)
    }else{
      y = lfun(a = a[ind],b = b[ind],S = S[ind,ind],nu = nu)
    }
    return(y)
  }
  R = S/(s%*%t(s))

  if(nu>0){
    a1s = sqrt(nu)*as
    b1s = sqrt(nu)*bs
    #print(nu)
    #y = pmvt(lower = a1s,upper = b1s,df = nu,corr = R)[1]
    log2y = pmvt.genz(lower = a1s,upper = b1s,nu = nu,sigma = R,uselog2 = TRUE)$Estimation
    y = gamma(nu/2)*2^log2y
  }else{
    if(nu == 0){
      if(sum(is.infinite(as - bs)) == n){ #at least one limit contains an Inf
        y = NaN
      }else{
        if(n==2){
          r = R[1,2]
          faux = function(w1){
            return((atan((bs[2]-r*w1)/sqrt((1-r^2)*(1+w1^2)))-atan((as[2]-r*w1)/sqrt((1-r^2)*(1+w1^2))))/sqrt(1+w1^2))
          }
          y = integrate(faux,as[1],bs[1],abs.tol = 1e-10,rel.tol = 1e-10)$value/pi
        }else{
          #if(n<=4){
          # ind1 = min(seq_len(n)[is.finite(as*bs)])  # find a finite pair of (a_i,b_i)
          # a2 = a[-ind1]
          # b2 = b[-ind1]
          # v1 = S[ind1,ind1]
          # c  = S[-ind1,ind1]/v1
          # R2 = S[-ind1,-ind1]-c%*%t(S[ind1,-ind1])
          # s2 = sqrt(diag(R2))
          # R2 = R2/(s2%*%t(s2))
          #
          # #integral
          # integd = function(x1){
          #   #n1 = length(x1)
          #   sx1 = sqrt(1+x1^2/v1)
          #   a2m = (a2-c*x1)/(s2*sx1)
          #   b2m = (b2-c*x1)/(s2*sx1)
          #   return(1/log2ratio(sx1,pmvt.genz(lower = c(a2m),upper = c(b2m),nu = 1,sigma = R2,uselog2 = TRUE)$Estimation))
          #   #return(pmvt.genz(lower = c(a2m),upper = c(b2m),nu = 1,sigma = R2)$Estimation/sx1)
          #     #pmvt(lower = c(a2m),upper = c(b2m),df = 1,corr = R2)[1]/sx1)
          # }
          #
          # integd2 = function(x2){
          #   return(sapply(x2,integd))
          # }
          # y = tryCatch(integrate(integd2,a[ind1],b[ind1],abs.tol = 1e-8,rel.tol = 1e-8)$value/sqrt(v1),
          #              error=function(e){
          #                seqq = seq(a[ind1],b[ind1],length.out = 1000);
          #                y = (b[ind1] - a[ind1])*mean(integd2(seqq))/sqrt(v1)})
          #}else{
          y = (lfun(a,b,S,-1e-1)+lfun(a,b,S,1e-1))/2
          #y = (lfun(a,b,S,-1)+lfun(a,b,S,1))/2   # When n>4, integration is not accurate, so use an approximation
          #}
        }
      }
    }else{  # dealing with nu<0
      y = 2/nu*lfun(a,b,S,nu+2)
      c = 1/(sqrt(pi)*nu)
      for(i in 1:n){
        a1 = as[i]
        b1 = bs[i]
        ma = S[-i,i]*a1/s[i]
        mb = S[-i,i]*b1/s[i]
        SS = S[-i,-i]-S[-i,i]%*%t(S[-i,i])/S[i,i]
        if(!is.infinite(a1)){
          y = y+c*a1/(1+a1^2)^((nu+1)/2)*lfun(a[-i]-ma,b[-i]-ma,(1+a1^2)*SS,nu+1)
        }
        if(!is.infinite(b1)){
          y = y-c*b1/(1+b1^2)^((nu+1)/2)*lfun(a[-i]-mb,b[-i]-mb,(1+b1^2)*SS,nu+1)
        }
      }
    }
  }
  return(y)
}
