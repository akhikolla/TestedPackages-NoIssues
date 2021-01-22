WenF<-function(pheRaw=NULL,genRaw=NULL,mapRaw1=NULL,yygg1=NULL,cov_en=NULL,
               WalkSpeed=NULL,CriLOD=NULL,dir=NULL){
  
  cl<-WalkSpeed;sLOD<-CriLOD;yygg<-NULL
  mx=NULL;phe=NULL;chr_name=NULL;v.map=NULL
  a.gen.orig=NULL;d.gen.orig=NULL;g1=NULL;X.ad.t4=NULL
  
  if(is.null(genRaw)==TRUE){
    showModal(modalDialog(title = "Warning", strong("Please input correct genotype dataset!"), easyClose = TRUE))
    
  }
  if(is.null(pheRaw)==TRUE){
    showModal(modalDialog(title = "Warning", strong("Please input correct phenotype dataset!"), easyClose = TRUE))
    
  }
  if(is.null(mapRaw1)==TRUE){
    showModal(modalDialog(title = "Warning", strong("Please input correct linkage map dataset!"), easyClose = TRUE))
    
  }
  if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(cl<0)){
    showModal(modalDialog(title = "Warning", strong("Please input Walk Speed: >0!"), easyClose = TRUE))
    
  }
  if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(cl>0)&&(sLOD<0)){
    showModal(modalDialog(title = "Warning", strong("Please input critical LOD score: >0!"), easyClose = TRUE))
    
  }
  
  if(is.null(yygg1)==FALSE){
    cov_en<-as.matrix(cov_en)
    yygg1<-as.matrix(yygg1)
    covname<-cov_en[2:nrow(cov_en),1]
    yygg1<-cbind(covname,yygg1)
    mapRaw10<-as.matrix(mapRaw1[-1,])
    chr_name<-unique(mapRaw10[,2])
    chr_secon<-as.matrix(mapRaw10[,2])
    mm<-numeric()
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_i<-which(chr_secon[]==chr_name[i])
      len<-matrix(length(chr_i),,1)
      mm<-rbind(mm,len)
      chr_name[i]<-i
    }
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_name<-as.numeric(chr_name)
      chr_pos1<-matrix(chr_name[i],mm[i],1)
      map_chr<-rbind(map_chr,chr_pos1)
    }
    map_marker<-matrix(mapRaw10[,1],,1)
    map_pos<-matrix(mapRaw10[,3],,1)
    mapRaw<-cbind(map_marker,map_chr,map_pos)
    nameMap<-matrix(mapRaw[,1],,1)
    nameGenrow<-matrix(genRaw[1,],1,)#individual's name in genotype 
    nameGencol<-matrix(genRaw[,1],,1)#marker's name in genotype
    namePhe<-as.matrix(pheRaw[,1],,1)
    nameCov<-matrix(yygg1[,1],,1)
    if(nameGenrow[2]=="1"){
      phee<-as.matrix(pheRaw[-1,-1])
      phe<-matrix(as.numeric(phee),nrow(phee),ncol(phee))
      mapname<-mapRaw
      genn<-genRaw[-1,-1]
      genn<-as.matrix(genn)
      mapnametwo<-matrix(as.numeric(mapname[,2:3]),nrow(mapname),2)
      mx<-cbind(mapnametwo,genn)
      chrRaw_name<-unique(mapRaw10[,2])
      newyygg<-as.matrix(yygg1)
      yyggChar<-as.matrix(newyygg[,-1])
      yygg<-matrix(as.numeric(yyggChar),nrow(yyggChar),ncol(yyggChar))
      mapname<-mapname
      chrRaw_name<-chrRaw_name
      chr_name<-chr_name
    }else{
      sameName_MG<-intersect(nameMap,nameGencol)
      sameName_PG<-intersect(namePhe,nameGenrow)
      locPhe<-match(sameName_PG,namePhe)
      locMap<-match(sameName_MG,nameMap)
      locGen_PG<-match(sameName_PG,nameGenrow)
      locGen_MG<-match(sameName_MG,nameGencol)
      locCov<-match(sameName_PG,nameCov)
      newyygg<-as.matrix(yygg1[locCov,])
      yyggChar<-as.matrix(newyygg[,-1])
      yygg<-matrix(as.numeric(yyggChar),nrow(yyggChar),ncol(yyggChar))
      newPhe<-as.matrix(pheRaw[locPhe,])
      newMap<-as.matrix(mapRaw[locMap,])
      newGenrow<-as.matrix(genRaw[,locGen_PG])
      newGen<-as.matrix(newGenrow[locGen_MG,])
      gen_two<-newMap[,2:3]
      genChar<-cbind(gen_two,newGen)
      if(ncol(newPhe)==2)
      {
        pheChar<-as.matrix(newPhe[,2:ncol(newPhe)])
      }else if(ncol(newPhe)>2){
        pheChar<-newPhe[,2:ncol(newPhe)]
      }
      chrRaw_name<-unique(mapRaw10[,2])
      chr_name<-chr_name
      mx<-as.matrix(genChar)
      phe<-matrix(as.numeric(pheChar),nrow(pheChar),ncol(pheChar))
      mapname<-newMap
    }
  } else{
    mapRaw10<-as.matrix(mapRaw1[-1,])
    chr_name<-unique(mapRaw10[,2])
    chr_secon<-as.matrix(mapRaw10[,2])
    mm<-numeric()
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_i<-which(chr_secon[]==chr_name[i])
      len<-matrix(length(chr_i),,1)
      mm<-rbind(mm,len)
      chr_name[i]<-i
    }
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_name<-as.numeric(chr_name)
      chr_pos1<-matrix(chr_name[i],mm[i],1)
      map_chr<-rbind(map_chr,chr_pos1)
    }
    map_marker<-matrix(mapRaw10[,1],,1)
    map_pos<-matrix(mapRaw10[,3],,1)
    mapRaw<-cbind(map_marker,map_chr,map_pos)
    nameMap<-matrix(mapRaw[,1],,1)
    nameGenrow<-matrix(genRaw[1,],1,)
    nameGencol<-matrix(genRaw[,1],,1)
    namePhe<-as.matrix(pheRaw[,1],,1)
    
    if(nameGenrow[2]=="1"){
      phee<-as.matrix(pheRaw[-1,-1])
      phe<-matrix(as.numeric(phee),nrow(phee),ncol(phee))
      mapname<-mapRaw
      genn<-genRaw[-1,-1]
      genn<-as.matrix(genn)
      mapnametwo<-matrix(as.numeric(mapname[,2:3]),nrow(mapname),2)
      mx<-cbind(mapnametwo,genn)
      chrRaw_name<-unique(mapRaw10[,2])
      chr_name<-chr_name
      mapname<-mapname
      chrRaw_name<-chrRaw_name
    }else{
      sameName_MG<-intersect(nameMap,nameGencol)
      sameName_PG<-intersect(namePhe,nameGenrow)
      locPhe<-match(sameName_PG,namePhe)
      locMap<-match(sameName_MG,nameMap)
      locGen_PG<-match(sameName_PG,nameGenrow)
      locGen_MG<-match(sameName_MG,nameGencol)
      newPhe<-as.matrix(pheRaw[locPhe,])
      newMap<-as.matrix(mapRaw[locMap,])
      newGenrow<-as.matrix(genRaw[,locGen_PG])
      newGen<-as.matrix(newGenrow[locGen_MG,])
      gen_two<-newMap[,2:3]
      genChar<-cbind(gen_two,newGen)
      if(ncol(newPhe)==2)
      {
        pheChar<-as.matrix(newPhe[,2:ncol(newPhe)])
      }else if(ncol(newPhe)>2){
        pheChar<-newPhe[,2:ncol(newPhe)]
      }
      chrRaw_name<-unique(mapRaw10[,2])
      chr_name<-chr_name
      mx<-as.matrix(genChar)
      phe<-matrix(as.numeric(pheChar),nrow(pheChar),ncol(pheChar))
      mapname<-newMap
    }
  }
  mapp1<-as.numeric(mapname[,2:3])
  mapp1<-matrix(mapp1,,2)
  chr<-length(unique(mapp1[,1]))
  for(i in 1:chr){
    pos1<-as.matrix(mapp1[which(mapp1[,1]==i),])
    delerow<-which(duplicated(pos1[,2]))
    if(length(delerow)!=0){
      break
    }else{
      mapname<-mapname
    }
  }
  if(length(delerow)!=0){
    showModal(modalDialog(title = "Warning", strong("Please check linkage maps (linkage groups) to make sure whether all the marker positions are different!"), easyClose = TRUE))
  }else{
    
    phe1<-as.matrix(phe[,1])
    mx<-as.matrix(mx)
    
    g1<-cbind(as.matrix(mapname[,1]),mx)
    blank<-as.matrix(c("phe","",""))
    p1<-t(rbind(blank,phe1))
    pgcombine<-rbind(p1,g1)
    
    # suppressWarnings(dir.create(path="temp",recursive = T))
    # dir1<-"temp/"
    
    write.table(pgcombine,file=paste(dir,"/listeria_rotY",".csv",sep=""),sep=",",row.names = F,col.names = F)
    
    Genotypes=c("A","H","B","D","C")
    Crosstype="f2"
    data.qtl<-read.cross("csvr",dir,"listeria_rotY.csv",genotypes=Genotypes,na.strings = "-",crosstype=Crosstype)
    data.qtl<-jittermap(data.qtl)
    
    bin<-cl  #user's options 1cM
    data.qtl.2<-calc.genoprob(data.qtl,step=bin,error.prob = 0.001)
    rownames(g1)<-NULL
    map.raw<-as.data.frame(g1[,1:3],stringsAsFactors = F)
    map.raw[,2:3]<-sapply(map.raw[,2:3],as.numeric)
    names(map.raw)<-c("marker","chr","pos")
    m0<-nrow(g1)
    
    gen.raw1<-cbind(mapRaw1,genRaw[,2:ncol(genRaw)])
    gen.raw<-gen.raw1[-1,]
    colnames(gen.raw)<-c("id","","",gen.raw1[1,4:ncol(gen.raw)])
    nchr<-max(as.numeric(mx[,1]))
    
    gen3<-calc.genoprob(data.qtl,step=0,error.prob = 0.001)
    marker.aa<-NULL
    marker.dd<-NULL
    for(rr in 1:nchr){
      mapgen<-gen3$geno[[rr]]$prob
      aam.gen<-mapgen[,,1]-mapgen[,,3]
      marker.aa<-cbind(marker.aa,aam.gen)
      ddm.gen<-mapgen[,,2]
      marker.dd<-cbind(marker.dd,ddm.gen)
    }
    a.gen.orig<-t(marker.aa)
    d.gen.orig<-t(marker.dd)
    
    aa.0<-NULL
    dd.0<-NULL
    v.map<-NULL
    nchr<-max(as.numeric(mx[,1]))
    for(ii in 1:nchr){
      map.gen<-data.qtl.2$geno[[ii]]$prob
      aa.gen<-map.gen[,,1]-map.gen[,,3]
      dd.gen<-map.gen[,,2]
      aa.0<-cbind(aa.0,aa.gen)
      dd.0<-cbind(dd.0,dd.gen)
      pos<-attr(map.gen,"map")
      pos.gen.1<-data.frame(pos)
      pos.gen.2<-data.frame(marker=row.names(pos.gen.1),chr=ii,pos=pos.gen.1,locus=1:dim(pos.gen.1)[1])
      id.insert<-grep(pattern = "loc",x=c(as.character(pos.gen.2$marker)))
      insertflag<-rep(1,dim(pos.gen.2)[1])
      insertflag[id.insert]<-0
      pos.gen.3<-data.frame(pos.gen.2,insertflag=insertflag,leftmarker=pos.gen.2$marker,rightmarker=pos.gen.2$marker)
      id.left<-findInterval(pos.gen.3[id.insert,]$pos,map.raw[map.raw$chr==ii,]$pos)
      namesleft<-map.raw[map.raw$chr==ii,]$marker[id.left]
      namesright<-map.raw[map.raw$chr==ii,]$marker[id.left+1]
      pos.gen.3[id.insert,]$leftmarker<-namesleft
      pos.gen.3[id.insert,]$rightmarker<-namesright
      pos.gen.3[id.insert,]$marker<-namesleft
      v.map<-rbind(v.map,pos.gen.3)
    }
    
    ##########################
    names.insert<-v.map
    ######################
    n<-ncol(g1)-3
    m.a<-dim(v.map)[1]
    #############################
    ad<-matrix(NA,nrow = (2*m.a),ncol=n)
    ad[seq(1,(2*m.a),by=2),]<-t(aa.0)
    ad[seq(2,(2*m.a),by=2),]<-t(dd.0)
    ad.gen.insert<-ad
    names.insert1<-names.insert[rep(1:nrow(names.insert),each=2),]
    names.insert2<-cbind(names.insert1,id.all=1:(2*m.a))
    X.ad.tran.data<-cbind(names.insert2,ad.gen.insert)
    X.ad.t4<-t(ad)
    
  }
  output<-list(yygg=yygg,mx=mx,phe=phe,chr_name=chr_name,v.map=v.map,gen.raw=gen.raw,a.gen.orig=a.gen.orig,d.gen.orig=d.gen.orig,n=n,
               names.insert2=names.insert2,X.ad.tran.data=X.ad.tran.data,X.ad.t4=X.ad.t4)
  return(output)
}



WenS<-function(flag=NULL,CriLOD=NULL,NUM=NULL,pheRaw=NULL,Likelihood=NULL,setseed=NULL,
               flagrqtl=NULL,yygg=NULL,mx=NULL,phe=NULL,chr_name=NULL,v.map=NULL,
               gen.raw=NULL,a.gen.orig=NULL,d.gen.orig=NULL,n=NULL,
               names.insert2=NULL,X.ad.tran.data=NULL,X.ad.t4=NULL,dir=NULL){
  
  sLOD<-CriLOD;result=NULL;mxmp=NULL;galaxyy1=NULL;res1a=NULL;res1d=NULL
  
  
  #########
  #-1,0,1 change a and d matrix
  change<-function(x){
    m<-length(x)
    y<-vector(length=2*m)
    for(j in 1:m){
      if(x[j]==1){y[(j-1)*2+1]<-1}
      if(x[j]==-1){y[(j-1)*2+1]<--1}
      if(x[j]==0){y[(j-1)*2+2]<-1}
    }
    return(y)
  }
  #####k kinship start################
  kinship.F2 <- function(gen){
    X1<-as.matrix(gen)
    kk<-X1%*%t(X1)
    cc<-mean(diag(kk))
    kk<-kk/cc
    return(kk)
  }
  ###################################################
  #K kinship end
  #############################################
  ###mixed.vars start
  #Mixednew can fit the models involving mutiple variace matrixs. 
  
  #input
  ##dataframe<-cbind(y,w)
  #y:n*1
  #w:n*c,include intercept,fixed effect
  #kinship<-list(K.a,K.d,...)
  #K.a,K.d,...,:n*n,matrix
  
  #output
  #RR$tau.kk:sigma.a2,sigma.d2,...,sigma.e2
  
  mixed.vars<-function(dataframe,kinship,optim.speed=TRUE){
    #data
    d<-dataframe
    y<-d[,1]
    x<-as.matrix(d[,-1,drop=FALSE])
    kk<-kinship
    indi<-nrow(d)
    r<-ncol(x)
    num.kk<-length(kk)
    
    cal.xx<-function(x,y,H=NULL){
      x<-as.matrix(x)
      y<-as.matrix(y)
      nr<-ncol(x)
      nc<-ncol(y)
      n0<-nrow(x)
      if (is.null(H)) H<-diag(n0)
      xx<-t(x)%*%H%*%y
      xx<-as.matrix(xx)
      return(xx)
    }
    
    
    ##likelihood function
    ##loglike
    loglike<-function(theta){
      i0<-num.kk+1
      V<-diag(indi)*as.numeric(theta[i0])
      for ( i in 1:num.kk){
        lambda<-as.numeric(theta[i])
        V<-V+kk[[i]]*lambda
      }
      
      v.inv<-solve(V,tol=1e-50)
      xx<-cal.xx(x=x,y=x,H=v.inv)
      xy<-cal.xx(x=x,y=y,H=v.inv)
      yy<-cal.xx(x=y,y=y,H=v.inv)
      
      d1<-unlist(determinant(V))
      d1<-d1[1]
      d2<-unlist(determinant(xx))
      d2<-d2[1]
      p<-yy-t(xy)%*%solve(xx)%*%xy
      
      
      logvalue<- -0.5*(d1+d2+p)
      return(-logvalue)
    }
    
    ## 
    ##parms involved in each kinship,sg2.1,sg2.2,...,se2. 
    if (optim.speed){
      theta<-rep(1,num.kk+1)
      parm<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=1e-10,upper=1e+10)
      tau.kk<-parm$par
      conv<-parm$convergence
      fn1<-parm$value
      hess<-parm$hessian
    }else{
      intervals<-c(1e-10,10,5e+02,1e+05)
      intials<-rep(1,num.kk+1)
      ni<-length(intervals)-1
      thetas<-list()
      for ( i in 1:ni){
        theta<-intials
        parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=intervals[1],upper=intervals[i+1])
        theta<-parms$par
        parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=1e-10,upper=1e+10)
        fi<-parms$par
        thetas[[i]]<-fi
      }
      fns<-sapply(1:ni,function(i){
        theta<-thetas[[i]]
        fn<-loglike(theta)
        return(fn)
      } )
      ii<-which.min(fns)
      tau.kk<-thetas[[ii]]
      fn1<-fns[ii]
    }
    
    ###  
    ###beta,var(beta),residual variance and blup estimation 
    i0<-num.kk+1
    V<-diag(indi)*as.numeric(tau.kk[i0])
    G<-matrix(0,indi,indi)
    for ( i in 1:num.kk){
      sg2<-as.numeric(tau.kk[i])
      V<-V+kk[[i]]*sg2
      G<-G+kk[[i]]*sg2
    }
    v.inv<-solve(V)
    xx<-cal.xx(x=x,y=x,H=v.inv)
    xy<-cal.xx(x=x,y=y,H=v.inv)
    yy<-cal.xx(x=y,y=y,H=v.inv)
    
    beta<-solve(xx,xy,tol=1e-50)
    v.beta<-solve(xx,tol=1e-50)      
    blup.sum<-G%*%v.inv%*%(y-x%*%beta)
    #Estimated weights
    yhat<-y-x%*%as.matrix(beta)
    yrandom<-NULL
    for (i in 1:num.kk){
      V0<-kk[[i]]*as.numeric(tau.kk[i])
      y0<-V0%*%v.inv%*%yhat
      yrandom<-cbind(yrandom,y0)
    }
    #
    xx<-t(yrandom)%*%yrandom
    xy<-t(yrandom)%*%yhat
    ww<-solve(xx,xy,tol=1e-50)
    ww<-as.vector(ww)
    #   
    
    RR<-list(beta=beta,v.beta=v.beta,random=blup.sum,se2=tau.kk[i0],tau.kk=tau.kk,weights=ww,lrt=fn1)
    ###
    return(RR)
  }
  #mixed.vars end
  ################################################
  #Y=w*alpha+Z*u+e start
  ################################################
  #likelihood
  emma.eigen.L <- function(Z,K,complete=TRUE) {
    if ( is.null(Z) ) {
      return(emma.eigen.L.wo.Z(K))
    }
    else {
      return(emma.eigen.L.w.Z(Z,K,complete))
    }
  }
  #likelihood
  emma.eigen.L.wo.Z <- function(K) {
    eig <- eigen(K,symmetric=TRUE)
    return(list(values=eig$values,vectors=eig$vectors))
  }
  #likelihood
  emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
    return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
  }
  
  #restricted likelihood
  emma.eigen.R <- function(Z,K,X,complete=TRUE) {
    if ( ncol(X) == 0 ) {
      return(emma.eigen.L(Z,K))
    }
    else if ( is.null(Z) ) {
      return(emma.eigen.R.wo.Z(K,X))
    }
    else {
      return(emma.eigen.R.w.Z(Z,K,X,complete))
    }
  }
  #restricted likelihood
  emma.eigen.R.wo.Z <- function(K, X) {
    n <- nrow(X)
    q <- ncol(X)
    S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
    eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
    stopifnot(!is.complex(eig$values))
    return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
  }
  #restricted likelihood
  emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
    if ( complete == FALSE ) {
      vids <-  colSums(Z) > 0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    n <- nrow(Z)
    t <- ncol(Z)
    q <- ncol(X)
    
    SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
    eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
    if ( is.complex(eig$values) ) {
      eig$values <- Re(eig$values)
      eig$vectors <- Re(eig$vectors)    
    }
    qr.X <- qr.Q(qr(X))
    return(list(values=eig$values[1:(t-q)],
                vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                             complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
  }
  
  emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
    n <- length(xi)
    delta <- exp(logdelta)
    return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) )  
  }
  
  emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
    delta <- exp(logdelta)
    return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
  }
  
  emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
    n <- length(xi)
    delta <- exp(logdelta)
    etasq <- etas*etas
    ldelta <- delta*lambda+1
    return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
  }
  
  emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
    delta <- exp(logdelta)
    etasq <- etas.1*etas.1
    ldelta <- delta*lambda+1
    return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
  }
  
  emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
    nq <- length(etas)
    delta <-  exp(logdelta)
    return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
  }
  
  emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
    tq <- length(etas.1)
    nq <- n - t + tq#(n-t):the number of eigvalue 1;tq:the number of components of etas.1
    delta <-  exp(logdelta)
    return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
  }
  
  emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
    nq <- length(etas)
    delta <- exp(logdelta)
    etasq <- etas*etas
    ldelta <- delta*lambda+1
    return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
  }
  
  emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
    t <- t1
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    etasq <- etas.1*etas.1
    ldelta <- delta*lambda+1
    return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
  }
  
  emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL)
  {
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    
    if ( det(crossprod(X,X)) == 0 ) {
      showModal(modalDialog(title = "Warning", strong("X is singular"), easyClose = TRUE))
      return (list(ML=0,delta=0,ve=0,vg=0))
    }
    
    if ( is.null(Z) ) {
      if ( is.null(eig.L) ) {
        eig.L <- emma.eigen.L.wo.Z(K)
      }
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.wo.Z(K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,n-q,m)    
      Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
      Xis.1<-matrix(eig.L$values,n,m)
      Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
      
      Etasq <- matrix(etas*etas,n-q,m)
      dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
        }
      }
      #    optdelta <- exp(optlogdelta)
    }
    else {
      if ( is.null(eig.L) ) {
        eig.L <- emma.eigen.L.w.Z(Z,K)
      }
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.w.Z(Z,K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      etas.1 <- etas[1:(t-q)]
      etas.2 <- etas[(t-q+1):(n-q)]
      etas.2.sq <- sum(etas.2*etas.2)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,t-q,m)
      Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
      
      Xis.1<-matrix(eig.L$values,t,m)
      Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
      Etasq <- matrix(etas.1*etas.1,t-q,m)
      
      dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
      #    optdelta <- exp(optlogdelta)
    }
    #handler of grids with NaN log
    optLL=replaceNaN(optLL)  #20160728
    
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    
    maxLL <- max(optLL)
    if ( is.null(Z) ) {
      maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
    }
    else {
      maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
    }
    maxvg <- maxve*maxdelta
    
    return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
    
  }
  
  
  emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                         esp=1e-10, eig.L = NULL, eig.R = NULL) {
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    
    if ( det(crossprod(X,X)) == 0 ) {
      showModal(modalDialog(title = "Warning", strong("X is singular"), easyClose = TRUE))
      return (list(REML=0,delta=0,ve=0,vg=0))
    }
    
    if ( is.null(Z) ) {
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.wo.Z(K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,n-q,m)
      Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
      Etasq <- matrix(etas*etas,n-q,m)
      dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
      #    optdelta <- exp(optlogdelta)
    }
    else {
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.w.Z(Z,K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      etas.1 <- etas[1:(t-q)]
      etas.2 <- etas[(t-q+1):(n-q)]
      etas.2.sq <- sum(etas.2*etas.2)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1 <- matrix(eig.R$values,t-q,m) 
      Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
      Etasq <- matrix(etas.1*etas.1,t-q,m)
      
      dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
      #    optdelta <- exp(optlogdelta)
    }  
    #handler of grids with NaN log
    optLL=replaceNaN(optLL)  #20160728
    
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    
    maxLL <- max(optLL)
    
    if ( is.null(Z) ) {
      maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
    }
    else {
      maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
    }
    maxvg <- maxve*maxdelta
    return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  }
  
  ################################################
  #Y=w*alpha+Z*u+e end
  
  #change to Y_c=W_c*alpha+X_c*beta+e_c,K=1 start
  ##################################################################################
  #change,k=1,Z=X_c,X=W_c:
  
  
  ################################################
  #likelihood 20190902
  FASTmrEMMA.delta.ML.LL.c<-function(logdelta,X,M,M.y,yMy,n){
    #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
    #n<-dim(M)[1]
    delta <-  exp(logdelta)
    ci<-as.numeric(crossprod(X))
    delta1<-as.numeric(t(X)%*%M%*%X)
    xMy<-as.numeric(crossprod(X,M.y))
    return(0.5*(n*((log(n/(2*pi))-log(as.numeric(yMy)-delta*(xMy)^2/(1+delta*delta1)))-1)-log(delta*ci+1)))
  }
  #dML
  
  FASTmrEMMA.delta.ML.dLL.c<-function(logdelta,X,M,M.y,yMy,n){
    #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
    #n<-dim(M)[1]
    delta <-  exp(logdelta)
    ci<-as.numeric(crossprod(X))
    delta1<-as.numeric(t(X)%*%M%*%X)
    xMy<-as.numeric(crossprod(X,M.y))
    return(-0.5*ci/(1+delta*ci)+0.5*n/((1+delta*delta1)*(as.numeric(yMy)*(1+delta*delta1)/(xMy^2)-delta)))
  }
  #restrict likelihood 20190902
  FASTmrEMMA.delta.REML.LL.c<-function(logdelta,X,M,M.y,yMy,v){
    #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
    #v<-n-1
    delta <-  exp(logdelta)
    #ci<-crossprod(X)
    delta1<-as.numeric(t(X)%*%M%*%X)
    xMy<-as.numeric(crossprod(X,M.y))
    return(0.5*(v*((log(v/(2*pi))-log(as.numeric(yMy)-delta*(xMy)^2/(1+delta*delta1)))-1)-log(delta*delta1+1)))
  }
  #dREML
  FASTmrEMMA.delta.REML.dLL.c<-function(logdelta,X,M,M.y,yMy,v){
    #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
    #n<-dim(M)[1]
    delta <-  exp(logdelta)
    #ci<-crossprod(X)
    delta1<-as.numeric(t(X)%*%M%*%X)
    xMy<-as.numeric(crossprod(X,M.y))
    return(-0.5*delta1/(1+delta*delta1)+0.5*v/((1+delta*delta1)*(as.numeric(yMy)*(1+delta*delta1)/(xMy^2)-delta)))
  }
  
  ######################################################################################
  #20190906
  FASTmrEMMA.MLE.c<-function(X,M,M.y,yMy,n, ngrids=100, llim=-10, ulim=10, esp=1e-10){
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    
    m <- length(logdelta)
    #delta <- exp(logdelta)
    
    dLL<-FASTmrEMMA.delta.ML.dLL.c(logdelta,X,M,M.y,yMy,n)
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL,FASTmrEMMA.delta.ML.LL.c(llim,X,M,M.y,yMy,n))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      #optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      optLL <- append(optLL, FASTmrEMMA.delta.ML.LL.c(ulim,X,M,M.y,yMy,n))
      
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        #r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        #r <- uniroot(FASTmrEMMA.delta.ML.dLL.c,lower = logdelta[i],upper = logdelta[i+1],X,M,M.y,yMy,n)
        r <- uniroot(FASTmrEMMA.delta.ML.dLL.c,c(logdelta[i],logdelta[i+1]),X=X,M=M,M.y=M.y,yMy=yMy,n=n)
        
        optlogdelta <- append(optlogdelta, r$root)
        #optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        optLL <- append(optLL,FASTmrEMMA.delta.ML.LL.c(r$root,X,M,M.y,yMy,n))
      }
    }
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    optLL=replaceNaN(optLL)
    maxLL <- max(optLL)
    
    xMy<-crossprod(X,M.y)
    xMx<-crossprod(X,(M%*%X))
    maxve <-(yMy-maxdelta*(xMy)^2/(1+maxdelta*xMx))/n
    #(sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
    
    maxvg <- maxve*maxdelta
    #alpha<-inv()
    
    return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,delta1=xMx,xMy=xMy))
  }
  
  FASTmrEMMA.REMLE.c<-function(X,M,M.y,yMy,v, ngrids=100, llim=-10, ulim=10, esp=1e-10){
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    
    m <- length(logdelta)
    #delta <- exp(logdelta)
    
    dLL<-FASTmrEMMA.delta.REML.dLL.c(logdelta,X,M,M.y,yMy,v)
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL,FASTmrEMMA.delta.REML.LL.c(llim,X,M,M.y,yMy,v))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      #optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      optLL <- append(optLL, FASTmrEMMA.delta.REML.LL.c(ulim,X,M,M.y,yMy,v))
      
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        #r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        r <- uniroot(FASTmrEMMA.delta.REML.dLL.c,lower = logdelta[i],upper = logdelta[i+1],X=X,M=M,M.y=M.y,yMy=yMy,v=v)
        
        optlogdelta <- append(optlogdelta, r$root)
        #optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        optLL <- append(optLL,FASTmrEMMA.delta.REML.LL.c(r$root,X,M,M.y,yMy,v))
      }
    }
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    optLL=replaceNaN(optLL)
    maxLL <- max(optLL)
    
    xMy<-crossprod(X,M.y)
    xMx<-crossprod(X,(M%*%X))
    maxve <-(yMy-maxdelta*(xMy)^2/(1+maxdelta*xMx))/v
    #(sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
    
    maxvg <- maxve*maxdelta
    #alpha<-inv()
    
    return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,delta1=xMx,xMy=xMy))
  }
  
  ###################################################################
  #################################################
  emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
    if( is.null(Z) ){
      return(emma.maineffects.B.Zo(K,deltahat.g))
    }
    else{
      return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
    }
  }
  
  #####
  emma.maineffects.B.Zo <-function(K,deltahat.g){
    t <- nrow(K)
    stopifnot(ncol(K) == t)
    
    B<-deltahat.g*K+diag(1,t)
    eig<-eigen(B,symmetric=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(mC=C,Q=Q,A=A))
  }
  
  emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    
    n <- nrow(Z)  
    B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
    eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(mC=C,Q=Q,A=A,complete=TRUE))
  }
  ##########################################################
  IRMMA.aK.dK.effects.B<-function(Z=NULL,aK,dK,deltahat.aK,deltahat.dK,complete=TRUE){
    if( is.null(Z)){
      return(IRMMA.aK.dK.effects.B.Zo(aK,dK,deltahat.aK,deltahat.dK))
    }
    else{
      return(IRMMA.aK.dK.effects.B.Z(Z,aK,dK,deltahat.aK,deltahat.dK,complete))
    }
  }
  ####
  IRMMA.aK.dK.effects.B.Zo<- function(aK,dK,deltahat.aK,deltahat.dK){
    n<- nrow(aK)
    stopifnot(nrow(dK)==n)
    
    stopifnot(ncol(dK)==n)
    stopifnot(ncol(dK)==n)
    
    B <- deltahat.aK*aK+deltahat.dK*dK+diag(1,n)
    eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(C.g=C,Q=Q,A=A))
  }
  
  ####
  IRMMA.aK.dK.effects.B.Z<-function(Z,aK,dK,deltahat.aK,deltahat.dK,complete){
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      aK <- aK[vids,vids]
      dK <- dK[vids,vids]
    }
    
    n <- nrow(Z)  
    t <- nrow(aK)
    
    stopifnot(ncol(aK)==t)
    stopifnot(nrow(dK)==t)
    
    stopifnot(ncol(dK)==t)
    
    B <- Z%*%(deltahat.aK*aK+deltahat.dK*dK)%*%t(Z)+diag(1,n)
    eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(C.g=C,Q=Q,A=A,complete=TRUE))
  }
  
  
  
  ################################################################################################
  emma.MLE0.c <- function(Y_c,W_c){
    
    n <- length(Y_c)
    
    stopifnot(nrow(W_c)==n)
    M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    etas<-crossprod(M_c,Y_c)
    LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
    return(list(ML=LL,M=M_c,n=n))
    
  }
  
  ######################################
  FASTmrEMMA.ML.LRT.c <- function(ys, xs, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
    #20190910
    #Z=C,X0=W=W0,xs=x:snp,n*p
    ys <- Z%*%ys   
    xs <- Z%*%xs
    X0 <- Z%*%X0
    
    ys<-as.matrix(ys)
    xs<-as.matrix(xs)
    X0<-as.matrix(X0)
    
    n <- nrow(ys)
    t <- ncol(xs)
    q<- if ( is.matrix(X0) ) ncol(X0) else 1
    v<-n-q
    
    MLE0<-emma.MLE0.c(ys,X0)
    
    ML1s <- vector(length=t)
    ML0s <- vector(length=t)
    vgs <- vector(length=t)
    ves <- vector(length=t)
    lambdas<-vector(length=t)
    bhats<-vector(length=t)
    #   
    d <- vector(length=t)
    stats <- vector(length=t)
    ps <- vector(length=t)
    
    #n<-199
    #M<-diag(1,n)-X0%*%ginv(crossprod(X0))%*%t(X0)
    M<-MLE0$M
    M.y<-M%*%ys
    yMy<-crossprod(ys,M.y)
    
    for (i in 1:t){
      
      # MLE1 <- emma.MLE.c (yv, x0v, K=1, xv, qr.X0,ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
      MLE1 <- FASTmrEMMA.MLE.c(X=xs[,i],M,M.y,yMy,n, ngrids=100, llim=-10, ulim=10, esp=1e-10)
      
      if(length(MLE1$vg)!=0){
        ML1s[i]<-MLE1$ML
        ML0s[i]<-MLE0$ML
        vgs[i]<-MLE1$vg
        ves[i]<-MLE1$ve
        lambdas[i]<-MLE1$delta
        
        ###################
        d[i] <- 1/(1+MLE1$delta*MLE1$delta1)
        #bhats[i]<-MLE1$lambda*MLE1$xMy/(1+MLE1$lambda*MLE1$delta1)
        bhats[i]<-MLE1$delta*MLE1$xMy/d[i]
        #to record me=sum(d) 
        
        
        stats[i]<- 2*(MLE1$ML-MLE0$ML)
        ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2#20160619
      }else{
        ps[i]<-1
      }
    }
    
    return(list(ID=1:t,ps=ps,bhats=bhats,lambdas=lambdas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
    
  } 
  
  emma.REMLE0.c <- function(Y_c,W_c){#20190831
    
    n <- length(Y_c)
    stopifnot(nrow(W_c)==n)
    
    t <-qr(W_c)$rank
    v <-n-t
    
    M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    etas<-crossprod(M_c,Y_c)
    
    LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
    return(list(REML=LL,M=M_c,v=v))
    
  }
  FASTmrEMMA.REML.LRT.c<- function(ys, xs, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
    
    
    ys <- Z%*%ys  
    xs <- Z%*%xs
    X0 <- Z%*%X0
    
    ys<-as.matrix(ys)
    xs<-as.matrix(xs)
    X0<-as.matrix(X0)
    
    n <- nrow(ys)
    t <- ncol(xs)
    q<- if ( is.matrix(X0) ) ncol(X0) else 1
    v<-n-q
    
    MLE0<-emma.REMLE0.c(ys,X0)
    
    ML1s <- vector(length=t)
    ML0s <- vector(length=t)
    vgs <- vector(length=t)
    ves <- vector(length=t)
    lambdas <- vector(length=t)
    bhats<-vector(length=t)
    
    d <- vector(length=t)
    
    stats <- vector(length=t)
    ps <- vector(length=t)
    
    M<-MLE0$M
    M.y<-M%*%ys
    yMy<-crossprod(ys,M.y)
    
    for (i in 1:t){
      
      #MLE1 <- emma.REMLE.c (ys, x0v, K=1, xv, qr.X0,ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)#20181112
      MLE1 <- FASTmrEMMA.REMLE.c(X=xs[,i],M,M.y,yMy,v, ngrids=100, llim=-10, ulim=10, esp=1e-10)
      
      if(is.na(MLE1$REML)==TRUE){
        ps[i]<-1
      }else{
        ML1s[i]<-MLE1$REML
        ML0s[i]<-MLE0$REML
        vgs[i]<-MLE1$vg
        ves[i]<-MLE1$ve
        lambdas[i] <- MLE1$delta
        ###################
        d[i] <- 1/(1+MLE1$delta*MLE1$delta1)
        #bhats[i]<-MLE1$lambda*MLE1$xMy/(1+MLE1$lambda*MLE1$delta1)
        bhats[i]<-MLE1$delta*MLE1$xMy/d[i]
        #to record me=sum(d) 
        
        stats[i]<- 2*(MLE1$REML-MLE0$REML)
        ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
      }
    }
    
    return(list(ID=1:t,ps=ps,bhats=bhats,lambdas=lambdas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vbs=vgs,ves=ves))
  }
  ######################################
  fixed.REMLE0.c <- function(Y_c,W_c){
    
    n <- length(Y_c)
    stopifnot(nrow(W_c)==n)
    M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    t <-qr(W_c)$rank
    v <-n-t
    etas<-crossprod(M_c,Y_c)
    
    LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
    return(list(REML=LL))
    
  }
  fixed.REML.LRT.c.sim<-function(ys=Y.orig[,i,drop=F], xs=X.ad.t4, Z=C2, X0=W.orig){
    
    ys <- Z%*%ys   
    xs <- Z%*%xs
    X0 <- Z%*%X0
    
    ys<-as.matrix(ys)
    xs<-as.matrix(xs)
    X0<-as.matrix(X0)
    
    
    n <- nrow(ys)
    m <- nrow(xs)
    t <- ncol(xs)
    stats <- vector(length=t)
    ps <- vector(length=t)
    
    
    ML1s <- vector(length=t)
    ML0s <- vector(length=t)
    
    
    REML0.c<-fixed.REMLE0.c(ys,X0)
    
    for(i in 1:t){
      #i<-2
      vids <- !is.na(xs[,i])
      xv <- xs[vids,i]
      
      yv <- ys[vids]
      x0v<-X0[vids,]
      
      xv.new<-cbind(x0v,xv)
      ##################3
      REML1.c<-fixed.REMLE0.c(yv,xv.new)
      
      xx.inv<-solve(t(xv.new)%*%xv.new)
      betahat<-xx.inv%*%t(xv.new)%*%ys
      sigma.e2<-(t(ys)%*%ys-t(ys)%*%xv.new%*%betahat)/(n-2)
      
      stats[i]<-betahat[2]^2/(sigma.e2*xx.inv[2,2])
      #stats[i]<- 2*(REML1.c$REML-REML0.c$REML)
      ps[i]<-pchisq(stats[i],1,lower.tail=F)#or pf(stats[i],1,n-2,lower.tail = F)
      
      
      ML1s[i]<-REML1.c$REML
      ML0s[i]<-REML0.c$REML
      
      
    }
    return(list(ID=1:t,ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats))
    
  } 
  ##########################
  fixed.ML.LRT.c.sim<-function(ys=Y.orig[,i,drop=F], xs=X.ad.t4, Z=C2, X0=W.orig){
    
    ys <- Z%*%ys   
    xs <- Z%*%xs
    X0 <- Z%*%X0
    
    ys<-as.matrix(ys)
    xs<-as.matrix(xs)
    X0<-as.matrix(X0)
    
    
    n <- nrow(ys)
    m <- nrow(xs)
    t <- ncol(xs)
    
    stats <- vector(length=t)
    ps <- vector(length=t)
    
    
    ML1s <- vector(length=t)
    ML0s <- vector(length=t)
    
    REML0.c<-emma.MLE0.c(ys,X0)#ys=Y_c,X0=W_c
    
    for(i in 1:t){
      vids <- !is.na(xs[,i])
      xv <- xs[vids,i]
      
      yv <- ys[vids]
      x0v<-X0[vids,]
      
      xv.new<-cbind(x0v,xv)
      
      REML1.c<-emma.MLE0.c(yv,xv.new)
      
      stats[i]<- 2*(REML1.c$ML-REML0.c$ML)
      ps[i]<-pchisq(stats[i],1,lower.tail=F)
      
      ML1s[i]<-REML1.c$ML
      ML0s[i]<-REML0.c$ML
      
      
    }
    return(list(ID=1:t,ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats))
    
  } 
  
  
  ######################################
  replaceNaN<-  function(LL) {
    #handler of grids with NaN log 
    index=(LL=="NaN")
    if(length(index)>0) theMin=min(LL[!index])
    if(length(index)<1) theMin="NaN"
    LL[index]=theMin
    return(LL)    
  }
  
  #########################################################################
  peak.id<-function(Lod.temp){
    m<-length(Lod.temp)
    optids<-vector(length=0)
    if(Lod.temp[1]>Lod.temp[2])   optids<-append(optids,1)
    
    for(j in 2:(m-1)){
      if ((Lod.temp[j-1]<Lod.temp[j]) & (Lod.temp[j]>Lod.temp[j+1])) {
        optids<-append(optids,j)
      }
    }
    if(Lod.temp[m]>Lod.temp[m-1])  optids<-append(optids,m)
    return(optids)
  }
  ##############################
  #Multinormal distribution density function
  multinormal<-function(y,mean,sigma)
  {
    pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
    return (pdf_value)
  }
  #LOD value test #library(MASS)
  ######################################3
  likelihood.a.d.F2<-function(xxn,xxx,yn,bbo,intercept)
    #xxn:fix matrix;xxx:gene matrix;yn:pheno matrix;bbo:gene effect from adalasso
  {
    nq<-ncol(xxx)
    ns<-nrow(yn)
    at1<-0
    ww1<-as.matrix(which(abs(bbo)>1e-5))
    
    ww.a<-ww1[ww1%%2==1]
    ww.d<-ww1[ww1%%2==0]
    
    ww.a.new<-c(ww.a,ww.a+1)
    ww.d.new<-c(ww.d,ww.d-1)
    
    ww1.new<-union(ww.a.new,ww.d.new)
    
    ww1.new<-ww1.new[order(ww1.new)]
    ww1.new<-as.matrix(ww1.new)
    
    at1<-dim(ww1.new)[1]
    lod<-matrix(rep(0,nq),nq,1)
    ps<-matrix(rep(1,nq),nq,1)
    
    ad<-if(at1>0.5) cbind(xxn,xxx[,ww1.new]) else xxn
    
    bb<-if (is.null(intercept)) ginv(crossprod(ad,ad))%*%crossprod(ad,yn) else c(intercept,bbo[ww1.new])
    
    vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns)
    ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
    
    sub<-1:ncol(ad)
    
    at2<-if(at1>1) seq(1,at1,by=2) else 1
    if(at1>0.5)
    {
      for(i in at2)
      {
        ij<-which((sub!=sub[i+ncol(xxn)])&(sub!=sub[i+ncol(xxn)+1]))
        
        ad1<-ad[,ij,drop=F]
        
        bb1<-if (is.null(intercept)) ginv(crossprod(ad1,ad1))%*%crossprod(ad1,yn) else c(intercept,bbo[ww1.new])[ij]
        
        vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns)
        ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
        lod[ww1.new[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
        
        ps[ww1.new[i]]<-pchisq(-2.0*(ll0-ll1),2,lower.tail = F)
      }
    }
    return (list(lod=lod,ps=ps))
  }
  
  adalasso<-function(X, y,k=10,use.Gram=TRUE,both=TRUE,intercept=TRUE){
    colnames(X)=1:ncol(X)
    n<-length(y)
    cv.adalasso<-NULL
    globalfit<-mylars(X,y,k=k,use.Gram=use.Gram,normalize=TRUE,intercept=intercept)
    coefficients.lasso=globalfit$coefficients
    intercept.lasso=globalfit$intercept
    cv.lasso<-globalfit$cv.lasso
    lambda<-globalfit$lambda
    lambda.lasso<-globalfit$lambda.opt
    coefficients.adalasso=NULL
    lambda.adalasso<-intercept.adalasso<-NULL
    if (use.Gram==TRUE){
      type="covariance"
    }
    if (use.Gram==FALSE){
      type="naive"
    }
    if (both==TRUE){ 
      # cross-validation for adaptive lasso
      #set.seed(11001)
      all.folds <- split(sample(1:n),rep(1:k,length=n))
      residmat <- matrix(0, length(lambda), k)
      
      for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain<-X[-omit,,drop=FALSE]
        ytrain<-y[-omit]
        Xtest<-X[omit,,drop=FALSE]
        ytest<-y[omit]
        my.lars<-mylars(Xtrain,ytrain,k=k,normalize=TRUE,use.Gram=use.Gram,intercept=intercept)
        coef.lasso<-my.lars$coefficients
        weights <- 1/abs(coef.lasso[ abs(coef.lasso)>0 ])
        #cat(paste("-- non-zero weights ",length(weights),"\n"))
        if (length(weights)==0){
          residmat[,i]<-mean((mean(ytrain)-ytest)^2)
        }
        if (length(weights)==1){
          residmat[,i]=mean((ytest -my.lars$intercept - Xtest%*%coef.lasso)^2)
        }
        if (length(weights)>1){
          XXtrain <- Xtrain[ , names(weights), drop=FALSE]
          XXtest<-Xtest[ , names(weights), drop=FALSE]
          XXtrain <- scale(XXtrain, center=FALSE, scale=weights)
          XXtest<-scale(XXtest, center=FALSE, scale=weights)
          #cat(paste("ncol of XXtrain: ",ncol(XXtrain),"\n"))
          fit<-glmnet(XXtrain,ytrain,type.gaussian=type,standardize=FALSE,intercept=intercept)
          pred<-predict(fit, newx=XXtest, type = "response",s = lambda)
          if (length(omit) == 1){
            pred <- matrix(pred, nrow = 1)
          }
          residmat[, i] <- apply((ytest - pred)^2, 2, mean)
        }
      }
      cv <- apply(residmat, 1, mean)
      cv.adalasso<-min(cv)
      weights <- 1/abs(coefficients.lasso[ abs(coefficients.lasso)>0 ])
      coefficients.adalasso<-rep(0,ncol(X))
      names(coefficients.adalasso)<-1:ncol(X)
      if (length(weights)>0){
        XX <- X[ , names(weights), drop=FALSE]
        if ( length(weights)==1 )  XX <- XX/weights        
        else  XX <- scale(XX, center=FALSE, scale=weights)
        if (length(weights)<=1){
          intercept.adalasso=intercept.lasso 
          coefficients.adalasso<-coefficients.lasso
          lambda.adalasso=0
        }
        else{
          fit<-glmnet(XX,y,type.gaussian=type,standardize=FALSE,intercept=intercept)
          lambda.adalasso<-lambda[which.min(cv)]
          coefficients=predict(fit,type="coefficients",s=lambda.adalasso)
          intercept.adalasso<-coefficients[1]
          coefficients.adalasso[names(weights)]<-coefficients[-1]/weights
        }
      }
    }
    return(list(cv.lasso=cv.lasso,lambda.lasso=lambda.lasso,cv.adalasso=cv.adalasso,lambda.adalasso=lambda.adalasso,intercept.lasso=intercept.lasso, intercept.adalasso=intercept.adalasso, coefficients.lasso=coefficients.lasso,coefficients.adalasso=coefficients.adalasso))
  }
  
  mylars<-function (X, y, k = 10,use.Gram=TRUE,normalize=TRUE,intercept=TRUE) 
  {
    x<-X
    n<-length(y)
    # set.seed(11001)
    all.folds <- split(sample(1:n),rep(1:k,length=n))
    
    if (use.Gram==TRUE){
      type="covariance"
    }
    if (use.Gram==FALSE){
      type="naive"
    }
    globalfit<-glmnet(x,y,family="gaussian",standardize=normalize,type.gaussian=type,intercept=intercept)
    lambda<-globalfit$lambda
    residmat <- matrix(0, length(lambda), k)
    for (i in seq(k)) {
      omit <- all.folds[[i]]
      fit <- glmnet(x[-omit, ,drop=FALSE], y[-omit],type.gaussian=type,standardize=normalize,family="gaussian",intercept=intercept)
      fit <- predict(fit, newx=x[omit, , drop = FALSE], type = "response", 
                     s = lambda)
      if (length(omit) == 1) 
        fit <- matrix(fit, nrow = 1)
      residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.lasso<-min(cv)
    cv.error <- sqrt(apply(residmat, 1, var)/k)
    lambda.opt<-lambda[which.min(cv)]
    coefficients=predict(globalfit,type="coefficients",s=lambda.opt)
    inter=coefficients[1]
    coefficients=coefficients[-1]
    names(coefficients)=1:ncol(X)
    object <- list(lambda=lambda,cv=cv,lambda.opt=lambda.opt,cv.lasso=cv.lasso,intercept=inter,coefficients=coefficients)
    invisible(object)
  }
 
   
    mxmp<-mx[,1:2]
    rownames(mxmp)<-NULL
    mxmp<-as.data.frame(mxmp,stringsAsFactors = F)
    mxmp[,1:2]<-sapply(mxmp[,1:2],as.numeric)
    
    m.a<-dim(v.map)[1]
    Y.orig<-as.matrix(phe[,NUM])
    Y.name<-pheRaw[1,NUM+1]
    ######################################
    
    W.orig<-matrix(1,n,1) 
    id.ind<-which(!is.na(Y.orig))
    Y.part.1<-Y.orig[id.ind,,drop=F]
    
    if(is.null(yygg)==FALSE){
      xenvir<-cbind(matrix(1,n,1),yygg) 
      xenvir<-xenvir[id.ind,]
      beta<-solve(t(xenvir)%*%xenvir)%*%t(xenvir)%*%Y.part.1
      Y.part.1<-Y.part.1-xenvir%*%beta+W.orig[id.ind,]
    }
    bb.adalasso<-NULL
    #kinship matrix#################
    a.K<-kinship.F2(t(a.gen.orig[,id.ind]))
    d.K<-kinship.F2(t(d.gen.orig[,id.ind]))
    kinship<-list(a.K,d.K)
    data.wy<-cbind(Y.part.1,rep(1,times=length(Y.part.1)))
    vars<-suppressWarnings(mixed.vars(data.wy,kinship,optim.speed=TRUE))
    delta.aK.init<-vars$tau.kk[1]/vars$tau.kk[3]
    delta.dK.init<-vars$tau.kk[2]/vars$tau.kk[3]
    
    remle.B<-IRMMA.aK.dK.effects.B(Z=NULL,a.K,d.K,delta.aK.init,delta.dK.init,complete=TRUE)
    C2<-remle.B$C.g
    if(flag==1){
      if (Likelihood=="REML"){
        REML.LRT.c2<-FASTmrEMMA.REML.LRT.c(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,],Z=C2, X0=W.orig[id.ind,,drop=F], ngrids=100, llim=-10, ulim=10, esp=1e-10)
      }else{
        REML.LRT.c2<-FASTmrEMMA.ML.LRT.c(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,],Z=C2, X0=W.orig[id.ind,,drop=F], ngrids=100, llim=-10, ulim=10, esp=1e-10)
      }
    }else{
      if(Likelihood=="REML"){
        REML.LRT.c2<-fixed.REML.LRT.c.sim(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,], Z=C2, X0=W.orig[id.ind,,drop=F])
      }else{
        REML.LRT.c2<-fixed.ML.LRT.c.sim(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,], Z=C2, X0=W.orig[id.ind,,drop=F])
      }
    }
    REML.LRT.c3<-data.frame(REML.LRT.c2)
    id.odd<-which(REML.LRT.c3$ID%%2==1)
    id.even<-which(REML.LRT.c3$ID%%2==0)
    LOD.a<--log10(REML.LRT.c3[id.odd,]$ps)
    LOD.d<--log10(REML.LRT.c3[id.even,]$ps)
    #########################3
    optids.a<-peak.id(LOD.a)
    optids.d<-peak.id(LOD.d)
    
    res1a<-cbind(v.map[,2:3],as.matrix(LOD.a))
    res1d<-cbind(v.map[,2:3],as.matrix(LOD.d))
    ###############################################
    ################################33
    a.eb.loc.1<-optids.a*2-1
    a.eb.loc.2<-c(a.eb.loc.1,a.eb.loc.1+1)
    d.eb.loc.1<-optids.d*2
    d.eb.loc.2<-c(d.eb.loc.1,d.eb.loc.1-1)
    ad.eb.loc.0<-union(a.eb.loc.2,d.eb.loc.2)
    ad.eb.loc.0<-ad.eb.loc.0[order(ad.eb.loc.0)]#xdata,a+d
    #########################################3
    xdata.ad<-X.ad.t4[id.ind,ad.eb.loc.0]
    X.ad.eb.names.0<-X.ad.tran.data[ad.eb.loc.0,1:8]
    
    suppressWarnings(RNGkind(sample.kind = "Rounding"))
    set.seed(setseed)  #user's options
    model1<-adalasso(X=xdata.ad,y=Y.part.1,k=10,use.Gram = T,both=T)
    adalasso.beta<-model1$coefficients.adalasso
    
    bb<-matrix(0,nrow=2*m.a,ncol=1)
    bb[ad.eb.loc.0,]<-adalasso.beta
    
    lrt.0<-likelihood.a.d.F2(xxn=W.orig[id.ind,,drop=F],xxx=xdata.ad,yn=Y.part.1,bbo=bb[ad.eb.loc.0,],intercept=model1$intercept.adalasso)
    id.lrt.0.a<-which(lrt.0$lod>=sLOD) #
    
    if(length(id.lrt.0.a)==0){
      showModal(modalDialog(title = "Warning", strong("No QTL were detected!"), easyClose = TRUE))
    }else if(length(id.lrt.0.a)>0){
      
      
      id.lrt.0.d<-id.lrt.0.a+1
      id.lrt.0.ad<-c(id.lrt.0.a,id.lrt.0.d)
      id.lrt.0.ad<-id.lrt.0.ad[order(id.lrt.0.ad)]
      
      opt.adalasso<-data.frame(X.ad.eb.names.0[id.lrt.0.a,],effect.a=bb[ad.eb.loc.0,][id.lrt.0.a],effect.d=bb[ad.eb.loc.0,][id.lrt.0.d],LOD=lrt.0$lod[id.lrt.0.a],ps=lrt.0$ps[id.lrt.0.a])
      
      ##########################
      if(flagrqtl=="TRUE"){
        
        # suppressWarnings(dir.create(path="temp",recursive = T))
        # dir1<-"temp/"
        
        b.qtl.1<-bb[ad.eb.loc.0,][id.lrt.0.ad]
        b.id.qtl<-which(b.qtl.1!=0)
        X.qtl<-X.ad.t4[id.ind,ad.eb.loc.0][,id.lrt.0.ad][,b.id.qtl]
        XY.qtl<-data.frame(Y.part.1,X.qtl)
        model.qtl.lm<-lm(Y.part.1~X.qtl,data = XY.qtl)
        
        b.qtl.2<-matrix(0,nrow = length(b.qtl.1),ncol=1)
        b.qtl.2[b.id.qtl]<-model.qtl.lm$coefficients[-1]
        
        ######################
        
        dif.1<-diff(opt.adalasso$locus)
        dif.2<-diff(opt.adalasso$pos)
        
        pos.options<-5 
        locus.options<-6 
        
        
        dif.id<-which(dif.1>0&dif.1<=locus.options&dif.2>0&dif.2<=pos.options)
        dif.id.ad.0<-c(dif.id,dif.id+1)
        
        dif.id.ad.0<-dif.id.ad.0[order(dif.id.ad.0)]
        qtl.id<-opt.adalasso[dif.id.ad.0,]
        dif.id.n<-length(dif.id.ad.0)
        dif.id.ad.1<-c()
        if(dif.id.n!=0){
          opt.qtl.all<-data.frame()
          for(ii in seq(1,dif.id.n,by=2)){
            if(((qtl.id[ii,]$effect.a==0&qtl.id[ii,]$effect.d!=0)&(qtl.id[ii+1,]$effect.a!=0&qtl.id[ii+1,]$effect.d==0))
               |((qtl.id[ii,]$effect.a!=0&qtl.id[ii,]$effect.d==0)&(qtl.id[ii+1,]$effect.a==0&qtl.id[ii+1,]$effect.d!=0))
            ){
              dif.id.ad.1<-c(dif.id.ad.1,dif.id.ad.0[c(ii,ii+1)])
              
              b.qtl.3<-b.qtl.2
              id.omit<-c(2*c(dif.id.ad.0[c(ii,ii+1)])-1,2*c(dif.id.ad.0[c(ii,ii+1)]))
              b.qtl.3[id.omit]<-0
              Y.part.2<-Y.part.1-X.ad.t4[id.ind,ad.eb.loc.0][,id.lrt.0.ad]%*%b.qtl.3
              Y.part.3<-matrix(NA,nrow=n,ncol = 1)
              Y.part.3[id.ind,]<-Y.part.2
              
              im.phe.rot<-t(cbind(Y.part.3,1:n))
              row.names(im.phe.rot)<-c(Y.name,"id")
              
              write.table(im.phe.rot,file=paste(dir,"/",NUM,"_phe_rot_",(ii+1)/2,".csv",sep=""),sep=",",row.names = T,col.names = F)
              ################
              marker.options<-5
              gen.qtl.id.left<-which(gen.raw[,1]==qtl.id[ii,]$leftmarker)
              gen.qtl.id.right<-which(gen.raw[,1]==qtl.id[ii+1,]$rightmarker)
              gen.qtl<-gen.raw[c((gen.qtl.id.left-marker.options):(gen.qtl.id.right+marker.options)),]
              
              im.gen.rot<-rbind(c("id","","",1:n),gen.qtl)
              write.table(im.gen.rot,file=paste(dir,"/",NUM,"_gen_rot_",(ii+1)/2,".csv",sep=""),sep=",",row.names = F,col.names = F)
              
              #####################################
              #rqtl
              data.qtl<-read.cross("csvsr",dir,paste(NUM,"_gen_rot_",(ii+1)/2,".csv",sep=""),paste(NUM,"_phe_rot_",(ii+1)/2,".csv",sep=""))
              data.qtl.1<-jittermap(data.qtl)
              data.qtl.2<-calc.genoprob(data.qtl.1,step=1,error.prob = 0.001)
              
              out.em<-suppressWarnings(scanone(data.qtl.2,method="em"))
              
              out1<-summary(out.em,threshold = 2.5)
              
              data.qtl.3<-sim.geno(data.qtl.1,step=1,n.draws = 1,error.prob =0.001 )
              qtl<-makeqtl(data.qtl.3,chr =out1$chr,pos=out1$pos )
              out2<-suppressWarnings(fitqtl(data.qtl.3,qtl=qtl,formula = y~Q1,dropone=F,get.ests = T))
              
              out3<-summary(out2)
              row.names(out3$ests)[2]
              
              opt.qtl<-data.frame(id.1=row.names(out1),id.2=row.names(out3$ests)[2],out1[,1:2],effect.a=out3$ests[2],effect.d=out3$ests[3],LOD=out1[,3],ps=out3$result.full[1,6],stringsAsFactors = F)
              
              if(opt.qtl$id.1%in%gen.raw[,1]){
                
                names.qtl.1<-names.insert2[which(names.insert2$chr==as.integer(as.character(opt.qtl$chr))),]
                
                names.left<-findInterval(opt.qtl$pos,names.qtl.1$pos)
                names.qtl.2<-names.qtl.1[names.left-1,]
                opt.qtl.lr<-data.frame(names.qtl.2,effect.a=opt.qtl$effect.a,effect.d=opt.qtl$effect.d,LOD=opt.qtl$LOD,ps=opt.qtl$ps,stringsAsFactors=F)
                
              }else{
                
                gen.qtl.1<-gen.qtl[which(gen.qtl[,2]==as.integer(as.character(opt.qtl$chr))),]
                left1<-findInterval(opt.qtl$pos,gen.qtl.1[,3])
                
                #################
                names.qtl.1<-names.insert2[which(names.insert2$chr==as.integer(as.character(opt.qtl$chr))),]
                names.left<-findInterval(opt.qtl$pos,names.qtl.1$pos)
                names.qtl.2<-names.qtl.1[names.left,]
                opt.qtl.lr<-data.frame(names.qtl.2[,1:2],locus=NA,pos=opt.qtl$pos,insertflag=0,leftmarker=names.qtl.2$leftmarker,rightmarker=gen.qtl.1[left1+1,1],id.all=names.qtl.2$id.all-1,effect.a=opt.qtl$effect.a,effect.d=opt.qtl$effect.d,LOD=opt.qtl$LOD,ps=opt.qtl$ps,stringsAsFactors=F)
                
              }
              opt.qtl.all<-rbind(opt.qtl.all,opt.qtl.lr)
            }
          }
          opt.im.all<-rbind(opt.adalasso[-dif.id.ad.1,],opt.qtl.all)
          opt.im.all.0<-opt.im.all[order(opt.im.all$chr,opt.im.all$pos),]
          
        }else{
          opt.im.all.0<-opt.adalasso
        }#if (length(dif.id.ad.0)!=0) end
        
      }else{
        opt.im.all.0<-opt.adalasso
      }
      ###################################
      X.r.id<-c(opt.im.all.0$id.all,opt.im.all.0$id.all+1)
      X.r.id<-X.r.id[order(X.r.id)]
      
      X.r<-X.ad.t4[id.ind,X.r.id]
      num.b<-length(X.r.id)
      b.r<-matrix(0,nrow =num.b)
      b.r[seq(1,num.b,by=2),]<-opt.im.all.0$effect.a
      b.r[seq(2,num.b,by=2),]<-opt.im.all.0$effect.d
      
      X.intc<-matrix(1,nrow=length(id.ind),ncol=1)
      intc<-model1$intercept.adalasso
      
      sigma.e2<-as.numeric(crossprod(Y.part.1-X.r%*%b.r-X.intc%*%intc,Y.part.1-X.r%*%b.r-X.intc%*%intc)/length(id.ind))
      
      r2.g<-as.numeric(sum(0.5*opt.im.all.0$effect.a^2+0.25*opt.im.all.0$effect.d^2))
      ge.all<-max(var(Y.part.1),(r2.g+sigma.e2))
      r2.p<-r2.g/ge.all
      r2.new<-(0.5*opt.im.all.0$effect.a^2+0.25*opt.im.all.0$effect.d^2)/ge.all
      
      Genei<-0.5*opt.im.all.0$effect.a^2+0.25*opt.im.all.0$effect.d^2
      vare<- c(round(sigma.e2,4),matrix("",length(Genei)-1,1))
      varp<- c(round(var(Y.part.1),4),matrix("",length(Genei)-1,1))
      
      #######################
      xdata.opt.adalasso<-data.frame(opt.im.all.0,r2=r2.new*100,PVE=r2.p*100,Genei,vare,varp)
      colnames(xdata.opt.adalasso)[c(13,14)]<-c("r^2(%)","PVE(%)")
      aa.adalasso<-data.frame(nametrait=rep(Y.name,times=dim(opt.im.all.0)[1]),xdata.opt.adalasso,stringsAsFactors=F)
      bb.adalasso<-aa.adalasso[,-c(2,5,6,9,13,15)]
      
      result<-bb.adalasso[,c(1,2,3,6,7,8,4,5,10,9,11,12)]
      rownames(result)<-NULL
      colnames(result)<-c("Trait","Chr","Position(cM)","Effect.a","Effect.d","LOD","Left_marker","right_marker","Var_Genet","r2(%)","Var_Error","Var_Phen(total)")
      
      galaxyy1<-as.matrix(result[,c(2,3,6)])
      
      # unlink(dir()[dir()=="temp"],recursive = T)
    }
    output<-list(result=result,mxmp=mxmp,galaxyy1=galaxyy1,res1a=res1a,res1d=res1d,chr_name=chr_name)
    return(output)
  }










