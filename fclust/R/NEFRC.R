
NEFRC = function(D, k, m, RS, startU,index,alpha,conv, maxit, seed = NULL)
{

  if (missing(D))
    stop("The distance matrix D must be given")
  if (is.null(D))
    stop("The distance matrix D is empty")
  D=data.matrix(D)
  n=nrow(D)
  p=ncol(D)
  distCheck(D,n,p)
  if (is.null(rownames(D)))
    rn=paste("Obj",1:n,sep=" ")
  else
    rn=rownames(D)
  if (is.null(colnames(D)))
    cn=paste("Var",1:p,sep=" ")
  else
    cn=colnames(D)
  if (any(is.na(D)))
    stop("The distance matrix D must not contain NA values")
  if (!is.numeric(D))
    stop("The distance matrix D is not a numeric data.frame or matrix")

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

          if(all(index!=c("PC","PE","MPC","SIL","SIL.F"))){
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

          if(!any(index==c("PC","PE","MPC","SIL","SIL.F"))){
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
  if (conv<=0)
  {
    cat("The convergence criterion conv must be a (small) value >0: the default value conv=1e-9 will be used ",fill=TRUE)
    conv=1e-9
  }
  if (!is.numeric(conv))
  {
    cat("The convergence criterion conv is not numeric: the default value conv=1e-9 will be used ",fill=TRUE)
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

  crit.f <- rep(NA,nk)
  crit = 0
  for(c in 1:nk)
  {
    if ((check!=1))
    {

      main.temp <- mainnefrc_U(D = D,U = startU,m = m,n = n,k = k[c],index = index,alpha = alpha,conv = conv,maxit = maxit)

    }else{

      main.temp <- mainnefrc(D = D,m = m,n = n,k = k[c],index = index,alpha = alpha,rs = RS,maxit = maxit,conv = conv)

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

  names(crit.f) =  paste(index," ","k=",k,sep="")
  k = main$k
  rownames(U.opt)=rn
  colnames(U.opt)=paste("Clus",1:k,sep=" ")
  names(value)=paste("Start",1:RS,sep=" ")
  names(it)=names(value)
  names(k)=c("Number of clusters")
  names(m)=c("Parameter of fuzziness")
  clus = cl.memb(U.opt)


  out = list(U=U.opt,
             H=NULL,
             F=NULL,
             clus=clus,
             medoid=NULL,
             value=value,
             criterion = crit.f,
             iter=it,
             k=k,
             m=m,
             ent=NULL,
             b=NULL,
             vp=NULL,
             delta=NULL,
             stand=NULL,
             D = D,
             call=match.call())

  class(out)=c("fclust")
  return(out)
}

