FKM.gk.noise <- function (X, k, m, vp, delta, RS, stand, startU, index,alpha,conv, maxit, seed = NULL)
  {
    if (missing(X))
      stop("The data set must be given")
    if (is.null(X))
      stop("The data set X is empty")
    X=data.matrix(X)
    n=nrow(X)
    p=ncol(X)

    if(p==1){
      stop("Dataset X must contain more than one dimension")
    }
    if (is.null(rownames(X)))
      rn=paste("Obj",1:n,sep=" ")
    else
      rn=rownames(X)
    if (is.null(colnames(X)))
      cn=paste("Var",1:p,sep=" ")
    else
      cn=colnames(X)
    if (any(is.na(X)))
      stop("The data set X must not contain NA values")
    if (!is.numeric(X))
      stop("The data set X is not a numeric data.frame or matrix")


    if ((missing(startU)) || (is.null(startU)))
    {
      check=1
      checkMiss = missing(k) || is.null(k) || !is.numeric(k)
      if (checkMiss)
      {
        k= 2:6
        cat("The default value k=2:6 has been set ",fill=TRUE)
      }
      nk <- length(k)
      if(nk != 1){

        if(!checkMiss){
          checkK = any(k <= 1) | any((k-as.integer(k)) != 0)
          if(checkK)
          {
            k  = 2:6
            nk <- length(k)
            cat("The number of clusters k must be greter than 1, integer and numeric: the default value k=2:6 will be used ",fill=TRUE)

          }else{

            k <- sort(unique(k))
            nk <- length(k)

          }
        }

        if (missing(index))
        {
          index = "SIL.F"
          cat("The default index SIL.F has been set ",fill=TRUE)
        }else{
          if(length(index) != 1)
          {
            index = "SIL.F"
            cat("The index must be a single value: SIL.F will be used ",fill=TRUE)
          }else{

            if(all(index!=c("PC","PE","MPC","SIL","SIL.F","XB"))){
              index = "SIL.F"
              cat("No match found for the index name: SIL.F will be used ", fill = TRUE)

            }
          }

        }
        if(index == "SIL.F"){
          if (missing(alpha))
          {
            alpha=1
            cat("The default value alpha=1 has been set for computing SIL.F ",fill=TRUE)
          }else{
            if (!is.numeric(alpha))
            {
              alpha=1
              cat("The weighting coefficient alpha is not numeric: the default value alpha=1 will be used for computing SIL.F ",fill=TRUE)
            }
            if (alpha<0)
            {
              alpha=1
              cat("The number of clusters k must be non negative: the value alpha=1 will be used for computing SIL.F ",fill=TRUE)
            }
          }
        }else{
          alpha = 1
        }
      }else{
        nk <- 1
        if (missing(index))
        {
          index = "SIL.F"
        }else{
          if(length(index) != 1)
          {
            index = "SIL.F"
            cat("The index must be a single value: SIL.F will be used ",fill=TRUE)
          }else{

            if(!any(index==c("PC","PE","MPC","SIL","SIL.F","XB"))){
              index = "SIL.F"
              cat("No match found for the index name: SIL.F will be used ", fill = TRUE)

            }
          }
        }

        if(index == "SIL.F"){
          if (missing(alpha))
          {
            alpha=1
          }else{
            if (!is.numeric(alpha))
            {
              alpha=1
              cat("The weighting coefficient alpha is not numeric: the default value alpha=1 will be used for computing SIL.F ",fill=TRUE)
            }
            if (alpha<0)
            {
              alpha=1
              cat("The number of clusters k must be non negative: the value alpha=1 will be used for computing SIL.F ",fill=TRUE)
            }
          }
        }else{
          alpha = 1
        }


        if (!is.numeric(k))
        {
          k=2
          cat("The number of clusters k is not numeric: the default value k=2 will be used ",fill=TRUE)
        }
        if ((k>ceiling(n/2)) || (k<2))
        {
          k=2
          cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: the default value k=2 will be used ",fill=TRUE)
        }
        if (k%%ceiling(k)>0)
        {
          k=ceiling(k)
          cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(nrow(X)/2)}: the value ceiling(k) will be used ",fill=TRUE)
        }
      }
    }else
    {
      startU=as.matrix(startU)
      ns=nrow(startU)
      k=ncol(startU)
      check=0
      nk = 1
      if (any(is.na(startU)))
      {
        k=2
        cat("The rational start must not contain NA values: the default value k=2 and a random start will be used ",fill=TRUE)
        check=1
      }
      if (!is.numeric(startU))
      {
        k=2
        cat("The rational start is not a numeric data.frame or matrix: the default value k=2 and a random start will be used ",fill=TRUE)
        check=1
      }
      if ((k>ceiling(n/2)) || (k<2))
      {
        k=2
        cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: the default value k=2 and a random start will be used ",fill=TRUE)
        check=1
      }
      if ((ns!=n) && (check=0))
      {
        cat("The number of rows of startU is different from that of X: k=ncol(startU) and a random start will be used ",fill=TRUE)
        check=1
      }
      if (any(apply(startU,1,sum)!=1))
      {
        startU=startU/apply(startU,1,sum)
        cat("The sums of the rows of startU must be equal to 1: the rows of startU will be normalized to unit row-wise sum ",fill=TRUE)
      }
      if (missing(index))
      {
        index = "SIL.F"
      }else{
        if(length(index) != 1)
        {
          index = "SIL.F"
          cat("The index must be a single value: SIL.F will be used ",fill=TRUE)
        }else{

          if(!any(index==c("PC","PE","MPC","SIL","SIL.F","XB"))){
            index = "SIL.F"
            cat("No match found for the index name: SIL.F will be used ", fill = TRUE)

          }
        }
      }

      if(index == "SIL.F"){
        if (missing(alpha))
        {
          alpha=1
        }else{
          if (!is.numeric(alpha))
          {
            alpha=1
            cat("The weighting coefficient alpha is not numeric: the default value alpha=1 will be used for computing SIL.F ",fill=TRUE)
          }
          if (alpha<0)
          {
            alpha=1
            cat("The number of clusters k must be non negative: the value alpha=1 will be used for computing SIL.F ",fill=TRUE)
          }
        }
      }else{
        alpha = 1
      }

    }
    if (missing(m))
    {
      m=2
    }
    if(length(m) != 1){
      m = 2
      cat("The parameter of fuzziness m must be a single value: the default value m=2 will be used",fill=TRUE)

    }
    if (!is.numeric(m))
    {
      m=2
      cat("The parameter of fuzziness m is not numeric: the default value m=2 will be used ",fill=TRUE)
    }
    if (m<=1)
    {
      m=2
      cat("The parameter of fuzziness m must be >1: the default value m=2 will be used ",fill=TRUE)
    }

    if (missing(vp))
    {if(nk > 1){
      vp = rep(1,max(k))
    }else{
      vp=rep(1,k)
    }
    }
    if (!is.numeric(vp))
    {
      if(nk > 1){
        cat("The volume parameter vp is not numeric: the default value vp=rep(1,max(k)) will be used ",fill=TRUE)
        vp = rep(1,max(k))
      }else{
        vp=rep(1,k)
        cat("The volume parameter vp is not numeric: the default value vp=rep(1,k) will be used ",fill=TRUE)
      }
    }
    if (!is.vector(vp))
    {
      if(nk > 1){
        cat("The volume parameter vp is not a vector: the default value vp=rep(1,max(k)) will be used ",fill=TRUE)
        vp = rep(1,max(k))
      }else{
        vp=rep(1,k)
        cat("The volume parameter vp is not a vector: the default value vp=rep(1,k) will be used ",fill=TRUE)
      }
    }
    if (min(vp)<=0)
    {
      if(nk > 1){
        cat("The volume parameter vp must be >0: the default value vp=rep(1,max(k)) will be used ",fill=TRUE)
        vp = rep(1,max(k))
      }else{
        vp=rep(1,k)
        cat("The volume parameter vp must be >0: the default value vp=rep(1,k) will be used ",fill=TRUE)
      }
    }
    if (length(vp)!=max(k))
    {
      if(nk > 1){
        cat("The number of elements of vp is different from max(k): the default value vp=rep(1,k) will be used ",fill=TRUE)
        vp = rep(1,max(k))
      }else{
        vp=rep(1,k)
        cat("The number of elements of vp is different from k: the default value vp=rep(1,k) will be used ",fill=TRUE)
      }
    }
    if (missing(RS))
    {
      RS=1
    }
    if (!is.numeric(RS))
    {
      cat("The number of starts RS is not numeric: the default value RS=1 will be used ",fill=TRUE)
      RS=1
    }
    if (RS<1)
    {
      cat("The number of starts RS must be an integer >=1: the default value RS=1 will be used ",fill=TRUE)
      RS=1
    }
    if (RS%%ceiling(RS)>0)
    {
      cat("The number of starts RS  must be an integer >=1: the value ceiling(RS) will be used ",fill=TRUE)
      RS=ceiling(RS)
    }
    if (missing(conv))
      conv=1e-9
    if (!is.numeric(conv))
    {
      cat("The convergence criterion conv is not numeric: the default value conv=1e-9 will be used ",fill=TRUE)
      conv=1e-9
    }
    if (conv<=0)
    {
      cat("The convergence criterion conv must be a (small) value >0: the default value conv=1e-9 will be used ",fill=TRUE)
      conv=1e-9
    }
    if (missing(maxit))
      maxit=1e+6
    if (!is.numeric(maxit))
    {
      cat("The maximum number of iterations maxit is not numeric: the default value maxit=1e+6 will be used ",fill=TRUE)
      maxit=1e+6
    }
    if (maxit<=0)
    {
      cat("The maximum number of iterations maxit must be an integer >0: the default value maxit=1e+6 will be used ",fill=TRUE)
      maxit=1e+6
    }
    if (maxit%%ceiling(maxit)>0)
    {
      cat("The maximum number of iterations maxit must be an integer >0: the value ceiling(maxit) will be used ",fill=TRUE)
      maxit=ceiling(maxit)
    }

    if(!is.null(seed))
    {
      if (!is.numeric(seed))
      {
        cat("The seed value is not numeric: set.seed(NULL) will be used ",fill=TRUE)
        set.seed(NULL)
      }else{
        set.seed(seed)
      }

    }

    Xraw=X
    rownames(Xraw)=rownames(X)
    colnames(Xraw)=colnames(X)
    if (missing(stand))
      stand=0
    if (!is.numeric(stand))
      stand=0
    if (stand==1)
      X=scale(X,center=TRUE,scale=TRUE)[,]
    checkDelta = FALSE
    if(missing(delta))
    {
      delta = rep(1,nk)
      checkDelta = TRUE

    }else{
      if(length(delta) != 1)
      {
        if(is.numeric(delta))
        {
          if(any(delta < 0)){
            cat("The noise distance delta must be non negative: the default value (see ?FKM.noise) will be used ",fill=TRUE)
            checkDelta = TRUE
            delta = rep(1,nk)
          }

          if(length(delta) != nk)
          {cat("The number of elements of delta is different from the number of elements of k: the default value (see ?FKM.noise) will be used" ,fill=TRUE)
            checkDelta = TRUE
            delta = rep(1,nk)

          }

          if (any(delta==0))
          {stop("When delta=0, the standard algorithm is applied: run the function FKM")}
        }else{
          checkDelta = TRUE
          delta = rep(1,nk)
          cat("The noise distance delta must is not numeric: the default value (see ?FKM.noise) will be used ",fill=TRUE)
        }

      }else{
        if (delta==0)
          stop("When delta=0, the standard algorithm is applied: run the function FKM")
        if (!is.numeric(delta))
        {
          cat("The noise distance delta must is not numeric: the default value (see ?FKM.noise) will be used ",fill=TRUE)
          checkDelta = TRUE
          delta = rep(1,nk)
        }
        if (delta<0)
        {
          cat("The noise distance delta must be non negative: the default value (see ?FKM.noise) will be used ",fill=TRUE)
          checkDelta = TRUE
          delta = rep(1,nk)
        }

        if(nk != 1 && !checkDelta)
        {
          delta = rep(delta,nk)
        }
      }
    }

    crit.f <- rep(NA,nk)
    crit = 0
    for(c in 1:nk)
    {
    if ((check!=1))
    {
      comp = k[c]
      if(checkDelta)
      {
        Hd=FKM.gk(X,k[c],m = m,RS=1,stand = stand,conv=1e-6,seed = seed)$H
        Dd = euclidean_distance(data = X,H = Hd,n = n,k = k[c])
        delta[c] <- sqrt(mean(Dd))
      }
      main.temp <- mainFKM_gk_noise_U(data = X, m = m,delta = delta[c], index = index,
                                      alpha = alpha,n = n, vp = vp[1:comp],p = p, k = k[c],
                                      U = startU, conv = conv, maxit = maxit)

      if(main.temp$warn == 1)
      {
        warning("At least one cluster covariance matrix seems to be singular. Change starting points or use FKM.gkb")
      }
    }else{
      comp = k[c]
      if(checkDelta)
      {
        Hd=FKM.gk(X,k[c],m = m,RS=1,stand = stand,conv=1e-6,seed = seed)$H
        Dd = euclidean_distance(data = X,H = Hd,n = n,k = k[c])
        delta[c] <- sqrt(mean(Dd))
      }
      main.temp <- mainFKM_gk_noise(data = X, m = m,delta = delta[c], index = index,
                                    alpha = alpha,n = n,vp = vp[1:comp] ,p = p, k = k[c],
                                    rs = RS, conv = conv, maxit = maxit)
      if(main.temp$warn == 1)
      {
        warning(paste("When k=", comp,", at least one cluster covariance matrix seems to be singular. Increase the number of starting points RS or use FKM.gkb",sep = ""))
      }
    }


    crit.temp = main.temp$index_max
    crit.f[c] = main.temp$index
    if(c == 1 | crit < crit.temp)
    {
      main = main.temp
      crit = crit.temp
    }
  }

value = as.vector(main$value)
it = as.vector(main$iter)

U.opt = main$U
H.opt = main$H
F.opt = main$F

names(crit.f) =  paste(index," ","k=",k,sep="")
names(delta) = paste("Noise distance k=",k)


k = main$k
vp = as.vector(main$vp)

rownames(H.opt)=paste("Clus",1:k)
colnames(H.opt)=cn
rownames(U.opt)=rn
colnames(U.opt)=rownames(H.opt)
dimnames(F.opt)[[1]]=cn
dimnames(F.opt)[[2]]=dimnames(F.opt)[[1]]
dimnames(F.opt)[[3]]=rownames(H.opt)
names(value)=paste("Start",1:RS,sep=" ")
names(it)=names(value)
names(k)=c("Number of clusters")
names(m)=c("Parameter of fuzziness")
names(vp)=rownames(H.opt)
if (stand!=1)
  stand=0
names(stand)=c("Standardization (1=Yes, 0=No)")
clus=cl.memb(U.opt)


out = list( U=U.opt,
            H=H.opt,
            F=F.opt,
            clus=clus,
            medoid=NULL,
            value=value,
            criterion = crit.f,
            iter=it,
            k=k,
            m=m,
            ent=NULL,
            b=NULL,
            vp=vp,
            delta=delta,
            stand=stand,
            Xca=X,
            X=Xraw,
            D=NULL,
            call=match.call())


    class(out)=c("fclust")
    return(out)
  }
