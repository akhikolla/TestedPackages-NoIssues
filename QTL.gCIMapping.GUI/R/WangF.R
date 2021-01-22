WangF<-function(pheRaw=NULL,genRaw=NULL,mapRaw1=NULL,yygg1=NULL,flagRIL=NULL,
                cov_en=NULL,Population=NULL,WalkSpeed=NULL,CriLOD=NULL){
  
  cl<-WalkSpeed;sLOD<-CriLOD;yygg<-NULL;chrRaw_name=NULL;
  mx=NULL;phe=NULL;chr_name=NULL;gen=NULL;mapname=NULL
  
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
  
  if(is.null(flagRIL)==TRUE){
    return(mapname)
  }else{
    
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
      showModal(modalDialog(title = "Warning", strong("Please check linkage maps (linkage groups) file to make sure whether all the marker positions are different!"), easyClose = TRUE))
    }else{
      mx<-as.matrix(mx)
      mx<-apply(mx,2,as.numeric)
      
      map<-mx[,1:2]
      geno<-t(mx[,3:(ncol(mx))])
      n_sam<-nrow(geno)
      gg1<-1;gg2<--1;gg0<-99
      
      
      mapinsert<-function(map,cl)
      {k<-0
      k1<-0
      mp<-numeric()
      for (ichr in 1:nrow(as.matrix(unique(map[,1]))))
      { 
        
        q1<-as.matrix(which(map[,1]==ichr))
        for (i in 2:(nrow(q1)))
        {rr<-map[q1[i],2]-map[q1[i-1],2]
        ll<-floor(rr/cl)
        q2<-rr-ll*cl
        if (q2>0){ll<-ll+1} 
        ss<-rr/ll
        k<-k+1
        for (j in 1:ll)
        {k1<-k1+1
        q3<-cbind((map[q1[i-1],2]+(j-1)*ss),ichr,(i-1),k,k1)
        mp<-rbind(mp,q3)
        }
        }
        k1<-k1+1
        q4<-cbind(map[q1[nrow(q1)],2],ichr,(nrow(q1)-1),k,k1)
        mp<-rbind(mp,q4)
      }
      return(mp)
      }
      
      mp<-mapinsert(map,cl)
      nq<-nrow(mp)
      mapp<-map
      mpp<-mp
      genoo<-geno
      
      
      markerall<-matrix(0,nrow(genoo),nrow(mpp))
      for (i in 1:nrow(as.matrix(unique(mapp[,1]))))
      {
        location<-as.matrix(which(mapp[,1]==i))
        pmap<-mapp[location,]
        pgeno<-genoo[,location]
        for (j in 1:(ncol(pgeno)-1))
        {
          ppmap<-pmap[j:(j+1),]
          ppgeno<-pgeno[,j:(j+1)]
          ppmap[,1]<-matrix(1,nrow(ppmap),1)
          map<-ppmap
          mp<-mapinsert(map,cl)
          ploc<-as.matrix(which(mpp[,2]==i & mpp[,1]>=map[1,2] & mpp[,1]<=map[2,2]))
          if (nrow(mp)>2)
          {for (ii in 1:n_sam)
          {geno<-t(as.matrix(ppgeno[ii,]))
          if ((geno[1]!=gg0) && (geno[2]!=gg0))
          {markerall[ii,ploc]<-QTL.gCIMapping.GUI::markerinsert(mp,geno,map,cl,gg1,gg2,gg0,flagRIL)
          }
          }
          }
          if (nrow(mp)==2)
          {w1<-which(ppgeno==gg1)
          w2<-which(ppgeno==gg2)
          ppgeno2<-ppgeno
          ppgeno2[w1]<-1
          ppgeno2[w2]<--1
          markerall[,ploc]<-ppgeno2
          }
        }
      }
      
      map<-mapp
      mp<-mpp
      geno<-genoo
      
      hm0<-matrix(0,nrow(geno),ncol(geno))
      for (i in 1:nrow(as.matrix(unique(map[,1]))))
      {
        
        hm<-as.matrix(which(map[,1]==i))
        for (j in 1:(nrow(hm)-1))
        {for (ii in 1:n_sam)
        {hmm<-as.matrix(cbind(geno[ii,hm[j]],geno[ii,hm[j+1]]))
        if (nrow(as.matrix(which(hmm==gg0)))==1)
        {if (as.matrix(which(hmm==gg0))==1)
        {hm0[ii,hm[j]]<-2
        hm0[ii,hm[j+1]]<-1
        }
          if (as.matrix(which(hmm==gg0))==2)
          {hm0[ii,hm[j]]<-1
          hm0[ii,hm[j+1]]<-2
          }
        }
        if (nrow(as.matrix(which(hmm==gg0)))==2)
        {hm0[ii,hm[j]]<-2
        hm0[ii,hm[j+1]]<-2
        }
        
        }
        }
      }
      
      for (i in 1:n_sam)
      {
        
        for (j in 1:nrow(as.matrix(unique(mapp[,1]))))
        {am<-as.matrix(which(mapp[,1]==j))
        pos0<-mapp[am,]
        hmm0<-t(as.matrix(hm0[i,am]))
        if (nrow(as.matrix(which(hmm0==1)))>0)
        {for (ii in 1:(nrow(as.matrix(which(hmm0==1)))+1))
        {amm<-as.matrix(which(hmm0==1))
        if (ii==1)
        {if (nrow(as.matrix(which(hmm0[1,1:amm[1]]==2)))>0)
        {geno<-t(as.matrix(genoo[i,am[1]:am[amm[1]]]))
        map<-pos0[1:amm[1],]
        loc<-as.matrix(which(mpp[,2]==map[1,1] & mpp[,1]>=map[1,2] & mpp[,1]<=map[nrow(map),2]))
        map[,1]<-map[,1]-map[1,1]+1
        mp<-mapinsert(map,cl)
        markerall[i,loc]<-QTL.gCIMapping.GUI::markerinsert(mp,geno,map,cl,gg1,gg2,gg0,flagRIL)
        }
        }
        if (ii>1 & ii<(nrow(as.matrix(which(hm0[i,am]==1)))+1))
        {if (nrow(as.matrix(which(hmm0[1,amm[ii-1]:amm[ii]]==2)))>0)
        {geno<-t(as.matrix(genoo[i,am[amm[ii-1]]:am[amm[ii]]]))
        map<-pos0[amm[ii-1]:amm[ii],]
        loc<-as.matrix(which(mpp[,2]==map[1,1] & mpp[,1]>=map[1,2] & mpp[,1]<=map[nrow(map),2]))
        map[,1]<-map[,1]-map[1,1]+1
        mp<-mapinsert(map,cl)
        markerall[i,loc]<-QTL.gCIMapping.GUI::markerinsert(mp,geno,map,cl,gg1,gg2,gg0,flagRIL)
        }
        }
        if (ii==(nrow(as.matrix(which(hm0[i,am]==1)))+1))
        {if (nrow(as.matrix(which(hmm0[1,amm[nrow(amm)]:ncol(hmm0)]==2)))>0)
        {geno<-t(as.matrix(genoo[i,am[amm[nrow(amm)]]:am[nrow(am)]]))
        map<-pos0[amm[nrow(amm)]:nrow(pos0),]
        loc<-as.matrix(which(mpp[,2]==map[1,1] & mpp[,1]>=map[1,2] & mpp[,1]<=map[nrow(map),2]))
        map[,1]<-map[,1]-map[1,1]+1
        mp<-mapinsert(map,cl)
        markerall[i,loc]<-QTL.gCIMapping.GUI::markerinsert(mp,geno,map,cl,gg1,gg2,gg0,flagRIL)
        }
        }
        }
        }
        }
      }
      map<-mapp
      mp<-mpp
      geno<-genoo
      gen<-cbind(mpp[,2],mpp[,1],t(markerall))
      
    }
    
    output<-list(yygg=yygg,mx=mx,phe=phe,chrRaw_name=chrRaw_name,chr_name=chr_name,gen=gen,mapname=mapname)
    return(output)
  }  
  
} 

#######################################################
WangS<-function(flag=NULL,CriLOD=NULL,NUM=NULL,pheRaw=NULL,chrRaw_name=NULL,
                yygg=NULL,mx=NULL,phe=NULL,chr_name=NULL,gen=NULL,mapname=NULL,CLO=NULL){
  
  sLOD<-CriLOD;result=NULL;mxmp=NULL;galaxyy1=NULL;res11=NULL;
  
  fix<-function(x,gen,y,kk){
    
    loglike<-function(theta){
      lambda<-exp(theta)
      logdt<-sum(log(lambda*delta+1))
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
      }
      if(abs(min(eigen(xx)$values))<1e-6)
        loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx+diag(ncol(xx))*0.01)%*%yx)-0.5*log(det(xx+diag(ncol(xx))*0.01))
      else
        loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)
    }
    
    fixed<-function(lambda){
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
      }
      if(abs(min(eigen(xx)$values))<1e-6)
        beta<-solve(xx+diag(ncol(xx))*0.01,yx)
      else
        beta<-solve(xx,yx)
      if(abs(min(eigen(xx)$values))<1e-6)
        sigma2<-(yy-t(yx)%*%solve(xx+diag(ncol(xx))*0.01)%*%yx)/(n-q)
      else
        sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      sigma2<-drop(sigma2)
      if(abs(min(eigen(xx)$values))<1e-6)
        vertue<-solve(xx+diag(ncol(xx))*0.01)
      else
        vertue<-solve(xx)
      var<-diag(vertue*sigma2)
      stderr<-sqrt(var)
      return(c(beta,stderr,sigma2))
    }
    qq<-eigen(as.matrix(kk))
    delta<-qq[[1]]
    uu<-qq[[2]]
    qx<-ncol(x)
    n<-length(y)
    yu<-t(uu)%*%y
    tempx<-x
    
    cl.cores <- detectCores()
    if(cl.cores<=2){
      cl.cores<-1
    }else if(cl.cores>2){
      if(cl.cores>10){
        cl.cores<-10
      }else {
        cl.cores <- detectCores()-1
      }
    }
    cll <- makeCluster(cl.cores)
    registerDoParallel(cll)
    
    i<-numeric()
    parmm<-foreach(i=1:nrow(gen),.combine=rbind)%dopar%{
      x<-tempx
      z<-gen[i,3:(ncol(gen))]
      qz<-ncol(z)
      x<-cbind(x,z)
      q<-ncol(x)
      xu<-t(uu)%*%x
      theta<-0
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
      lambda<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      lrt<-2*(fn0-fn1)
      hess<-parm$hessian
      parmfix<-fixed(lambda)
      beta<-parmfix[1:q]
      stderr<-parmfix[(q+1):(2*q)]
      sigma2<-parmfix[2*q+1]
      poly.lod<-lrt/4.61
      poly.p<-pchisq(lrt,1,lower.tail = F)
      sigma2g<-lambda*sigma2
      g<-beta[-c(1:qx)]
      g.err<-stderr[-c(1:qx)]
      b<-beta[c(1:qx)]
      b.err<-stderr[c(1:qx)]
      wald<-g^2/g.err^2
      p<-pchisq(wald,1,lower.tail = F)
      parmm<-(c(b[1],sigma2,lambda,sigma2g,poly.lod,poly.p,g,g.err,wald,p))
    }
    stopCluster(cll)
    return(parmm)
  }
  
  
  random<-function(fx,gen,phe,kk,CLO)
  {
    mixed<-function(x,y,kk){
      
      loglike<-function(theta){
        lambda<-exp(theta)
        logdt<-sum(log(lambda*delta+1))
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,q,1)
        xx<-matrix(0,q,q)
        for(i in 1:q){
          yx[i]<-sum(yu*h*xu[,i])
          for(j in 1:q){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])
          }
        }
        loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
        return(-loglike)
      }
      
      fixed<-function(lambda){
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,q,1)
        xx<-matrix(0,q,q)
        for(i in 1:q){
          yx[i]<-sum(yu*h*xu[,i])
          for(j in 1:q){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])
          }
        }
        beta<-solve(xx,yx)
        sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
        sigma2<-drop(sigma2)
        var<-diag(solve(xx)*sigma2)
        stderr<-sqrt(var)
        return(c(beta,stderr,sigma2))
      }
      
      qq<-eigen(kk)
      delta<-qq[[1]]
      uu<-qq[[2]]
      q<-ncol(x)
      n<-ncol(kk)
      vp<-var(y)
      yu<-t(uu)%*%y
      xu<-t(uu)%*%x
      theta<-0
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
      lambda<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      lrt<-2*(fn0-fn1)
      hess<-parm$hessian
      parmfix<-fixed(lambda)
      beta<-parmfix[1:q]
      stderr<-parmfix[(q+1):(2*q)]
      sigma2<-parmfix[2*q+1]
      lod<-lrt/4.61
      p_value<-pchisq(lrt,1,lower.tail = F)
      sigma2g<-lambda*sigma2
      goodness<-(vp-sigma2)/vp
      par<-data.frame(lrt,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
      return(par)
    }
    
    
    loglike<-function(theta){
      xi<-exp(theta)
      tmp0<-zz*xi+1
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      logdt2<-log(det(tmp0))
      loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
      return(-loglike)
    }
    
    fixed<-function(xi){
      tmp0<-zz*xi+diag(1)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      zHy<-zy-zz%*%tmp%*%zy
      zHx<-zx-zx%*%tmp%*%zz
      zHz<-zz-zz%*%tmp%*%zz
      beta<-solve(xHx,yHx)
      tmp2<-solve(xHx)
      sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
      gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
      var<-abs((xi*diag(1)-xi*zHz*xi)*as.numeric(sigma2))
      stderr<-sqrt(diag(var))
      result<-list(gamma,stderr,beta,sigma2)
      return(result)
    }
    name<-gen[,1:2]
    gen<-gen[,3:(ncol(gen))]
    gen<-t(gen)
    n<-nrow(gen)
    m<-ncol(gen)
    x<-fx
    
    s<-ncol(x)
    kk<-as.matrix(kk)
    qq<-eigen(kk)
    delta<-qq[[1]]
    uu<-qq[[2]]
    xu<-t(uu)%*%x
    y<-as.matrix(phe)
    parm<-mixed(x=x,y=y,kk=kk)
    lambda<-parm$lambda[1]
    h<-1/(delta*lambda+1)
    yu<-t(uu)%*%y
    xx<-matrix(0,s,s)
    for(i in 1:s){
      for(j in 1:s){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])
      }
    }
    yy<-sum(yu*h*yu)
    yx<-matrix(0,s,1)
    for(i in 1:s){
      yx[i]<-sum(yu*h*xu[,i])
    }
    if(is.null(CLO)==TRUE){
      cl.cores <- detectCores()
      if(cl.cores<=2){
        cl.cores<-1
      }else if(cl.cores>2){
        if(cl.cores>10){
          cl.cores<-10
        }else {
          cl.cores <- detectCores()-1
        }
      }
      cll <- makeCluster(cl.cores)
      registerDoParallel(cll)
      
      k<-numeric()
      parms<-foreach(k=1:m,.combine=rbind)%dopar%{
        z<-as.matrix(gen[,k])
        zu<-t(uu)%*%z
        zy<-as.matrix(sum(yu*h*zu))
        zz<-as.matrix(sum(zu*h*zu))
        zx<-matrix(0,s,1)
        for(i in 1:s){
          zx[i]<-sum(xu[,i]*h*zu)
        }
        theta<-c(0)
        par<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
        xi<-exp(par$par)
        conv<-par$convergence
        fn1<-par$value
        hess<-par$hessian
        parmfix<-fixed(xi)
        gamma<-parmfix[[1]]
        stderr<-parmfix[[2]]
        beta<-parmfix[[3]]
        sigma2<-parmfix[[4]]
        lambda<-xi
        sigma2g<-lambda*sigma2
        fn0<-loglike(-Inf)
        lrt<-2*(fn0-fn1)
        p_lrt<-pchisq(lrt,1,lower.tail = F)
        wald<-(gamma/stderr)^2
        p_wald<-pchisq(wald,1,lower.tail = F)
        parm0<-c(k,name[k,1],name[k,2],beta[1],sigma2,sigma2g,gamma,stderr,wald,p_wald)
      }
      stopCluster(cll)
    }else{
      
      qq<-numeric()
      for(k in 1:m){
        z<-as.matrix(gen[,k])
        zu<-t(uu)%*%z
        zy<-as.matrix(sum(yu*h*zu))
        zz<-as.matrix(sum(zu*h*zu))
        zx<-matrix(0,s,1)
        for(i in 1:s){
          zx[i]<-sum(xu[,i]*h*zu)
        }
        theta<-c(0)
        par<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
        xi<-exp(par$par)
        conv<-par$convergence
        fn1<-par$value
        hess<-par$hessian
        parmfix<-fixed(xi)
        gamma<-parmfix[[1]]
        stderr<-parmfix[[2]]
        beta<-parmfix[[3]]
        sigma2<-parmfix[[4]]
        lambda<-xi
        sigma2g<-lambda*sigma2
        fn0<-loglike(-Inf)
        lrt<-2*(fn0-fn1)
        p_lrt<-pchisq(lrt,1,lower.tail = F)
        wald<-(gamma/stderr)^2
        p_wald<-pchisq(wald,1,lower.tail = F)
        parm0<-c(k,name[k,1],name[k,2],beta[1],sigma2,sigma2g,gamma,stderr,wald,p_wald)
        qq<-rbind(qq,parm0)
      }
      parms<-qq
    }
    return(parms)
    
  }
  
  
  multinormal<-function(y,mean,sigma)
  {
    pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
    return (pdf_value)
  }
  #LOD value test
  likelihood<-function(xxn,xxx,yn)
  {
    nq<-ncol(xxx)
    ns<-nrow(yn)
    at1<-0
    ww1<-1:ncol(xxx)
    ww1<-as.matrix(ww1)
    at1<-dim(ww1)[1]
    lod<-matrix(rep(0,nq),nq,1)
    if(at1>0.5)
      ad<-cbind(xxn,xxx[,ww1])
    else
      ad<-xxn
    #if(abs(det(crossprod(ad,ad)))<1e-6)
    if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
      bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
    else
      bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
    vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
    ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
    
    sub<-1:ncol(ad);
    if(at1>0.5)
    {
      for(i in 1:at1)
      {
        ij<-which(sub!=sub[(i+ncol(xxn))])
        ad1<-ad[,ij]
        #if(abs(det(crossprod(ad1,ad1)))<1e-6)
        if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
          bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
        else
          bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn)
        vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
        ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
        lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
      }
    }
    return (lod)
  }
  
  #2010 EM_Bayes
  ebayes_EM<-function(x,z,y)
  {
    n<-nrow(z);k<-ncol(z)
    if(abs(min(eigen(crossprod(x,x))$values))<1e-6)
      b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
    else
      b<-solve(crossprod(x,x))%*%crossprod(x,y)
    v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
    u<-matrix(rep(0,k),k,1)
    v<-matrix(rep(0,k),k,1)
    s<-matrix(rep(0,k),k,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      s[i]<-((crossprod(zz,zz))^(-1))*v0
      u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
      v[i]<-u[i]^2+s[i]
    }
    vv<-matrix(rep(0,n*n),n,n);
    for(i in 1:k)
    {
      zz<-z[,i]
      vv=vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    iter<-0;err<-1000;iter_max<-100;err_max<-1e-8
    tau<-0;omega<-0
    while((iter<iter_max)&&(err>err_max))
    {
      iter<-iter+1
      v01<-v0
      v1<-v
      b1<-b
      vi<-solve(vv)
      xtv<-crossprod(x,vi)
      if(ncol(x)==1)
      {
        b<-((xtv%*%x)^(-1))*(xtv%*%y)
      }else
      {
        if(abs(min(eigen(xtv%*%x)$values))<1e-6){
          b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
        }
        else{
          b<-solve(xtv%*%x)%*%(xtv%*%y)
        }
      }
      r<-y-x%*%b
      ss<-matrix(rep(0,n),n,1)
      for(i in 1:k)
      {
        zz<-z[,i]
        zztvi<-crossprod(zz,vi)
        u[i]<-v[i]*zztvi%*%r
        s[i]<-v[i]*(1-zztvi%*%zz*v[i])
        v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
        ss<-ss+zz*u[i]
      }
      v0<-as.numeric(crossprod(r,(r-ss))/n)
      vv<-matrix(rep(0,n*n),n,n)
      for(i in 1:k)
      {
        zz<-z[,i]
        vv<-vv+tcrossprod(zz,zz)*v[i]
      }
      vv<-vv+diag(n)*v0
      err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
      beta<-t(b)
      sigma2<-v0
    }
    wang<-matrix(rep(0,k),k,1)
    for (i in 1:k){
      stderr<-sqrt(s[i]+1e-20)
      t<-abs(u[i])/stderr
      f<-t*t
      p<-pchisq(f,1,lower.tail = F)
      wang[i]<-p
    }
    return (wang)
  }
  
  phe<-as.matrix(phe[,NUM])
  
  if((is.null(pheRaw)==TRUE)&(is.null(chrRaw_name)==TRUE)&(is.null(yygg)==TRUE)&
     (is.null(mx)==TRUE)&(is.null(chr_name)==TRUE)&(is.null(gen)==TRUE)&(is.null(mapname)==TRUE)){
    return(phe)  
  }else{
    
    deletRow<-which(is.na(phe)==TRUE)
    gentwo<-mx[,1:2]
    t_gen<-mx[,3:ncol(mx)]
    
    if(length(deletRow)>0){
      phe<-as.matrix(phe[-deletRow,])
      t_gen1<-t_gen[,-deletRow]
      if(is.null(yygg)==FALSE){
        yygg<-yygg[-deletRow,]
      }
      phe<-phe
      mx<-cbind(gentwo,t_gen1)
      gen<-gen[,-(deletRow+2)]
    }else{
      mx<-mx
      phe<-phe
      if(is.null(yygg)==FALSE){
        yygg<-yygg
      }
      gen<-gen
    }
    
    mx<-as.matrix(mx)
    
    mxmp<-mx[,1:2]
    rownames(mxmp)<-NULL
    mxmp<-as.data.frame(mxmp,stringsAsFactors = F)
    mxmp[,1:2]<-sapply(mxmp[,1:2],as.numeric)
    
    
    map<-mx[,1:2]
    geno<-t(mx[,3:(ncol(mx))])
    n_sam<-nrow(geno)
    
    if(is.null(yygg)==FALSE){
      fx<-cbind(matrix(1,n_sam,1),yygg)
    }
    if(is.null(yygg)==TRUE){
      fx<-matrix(1,n_sam,1)
    }
    
    ori<-NULL
    for (j in 1:(nrow(map)))
    {
      ta<-as.matrix(which(gen[,1]==map[j,1]))
      tb<-gen[ta,2]
      cori<-ta[as.matrix(which(tb==map[j,2])),1]
      ori<-rbind(ori,cori)
    }
    corie<-matrix(0,nrow(gen),1)
    corie[ori,1]<-matrix(1,nrow(ori),1)
    gen<-cbind(gen[,1:2],corie,gen[,3:(ncol(gen))])
    wg<-gen
    iw<-as.matrix(which(gen[,3]==1))
    gk<-t(gen[iw,4:(ncol(gen))])
    m<-ncol(gk)
    n<-nrow(gk)
    kk<-matrix(0,n,n)
    for(k in 1:m){
      z<-as.matrix(gk[,k])
      kk<-kk+z%*%t(z)
    }
    cc<-mean(diag(kk))
    kk<-kk/cc
    gen<-cbind(gen[,1:2],gen[,4:(ncol(gen))])
    
    
    if (flag==1)
    {code<-random(fx=fx,gen=gen,phe=phe,kk=kk,CLO)
    tempcode<-code
    }
    if (flag==0)
    {code<-fix(x=fx,gen=gen,y=phe,kk=kk)
    tempcode<-code
    tempcode[,2:3] <- gen[,1:2]
    }
    res1<-tempcode
    
    res11<-res1[,c(2,3,10)]
    colnames(res11)<-NULL
    
    x0<-t(gen[,3:ncol(gen)])
    y<-phe
    bb<-code
    bb<-as.matrix(bb)
    aa<-numeric()
    for (i in 1:nrow(as.matrix(unique(gen[,1])))){
      
      mc<-which(gen[,1]==i)
      mc<-as.matrix(mc)
      for (j in 1:(nrow(mc)-2))
      {if (bb[mc[j+1],10]<bb[mc[j],10] & bb[mc[j+1],10]<bb[mc[j+2],10]){aa<-rbind(aa,mc[j+1])}
      }
      if (bb[mc[1],10]<bb[mc[2],10]){aa<-rbind(aa,mc[1])}
      if (bb[mc[nrow(mc)],10]<bb[mc[nrow(mc)-1],10]){aa<-rbind(aa,mc[nrow(mc)])}
    }
    
    xx<-x0[,aa]
    par<-ebayes_EM(fx,xx,y)
    selectpos<-which(par[,1]<=0.01)
    
    if(length(selectpos)==0){
      showModal(modalDialog(title = "Warning", strong("No QTL were detected!"), easyClose = TRUE))
    }else{
      cc<-which(par[,1]<=0.01)
      name<-as.matrix(aa[cc,1])
      xxx<-as.matrix(x0[,name])
      y<-as.matrix(y)
      lod<-likelihood(fx,xxx,y)
      dd<-as.matrix(which(lod[,1]>=sLOD))
      
      if(length(dd)==0){
        showModal(modalDialog(title = "Warning", strong("No QTL were detected!"), easyClose = TRUE))
      }else{
        
        if(length(dd)==1){
          na<-matrix(name[dd],1,)
          xxxm<-matrix(xxx[,dd],,1)
          wow<-cbind(fx,xxxm)
          bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
          ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          genna1<-matrix(gen[na,1],,1)
          genna2<-matrix(gen[na,2],,1)
          genna3<-as.matrix(ef)
          genna4<-as.matrix(lod[dd,])
          galaxy<-cbind(genna1,genna2,genna3,genna4)
        }else{
          na<-as.matrix(name[dd])
          wow<-cbind(fx,xxx[,dd])
          bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
          ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          galaxy<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
        }
        galaxy<-galaxy
        c1<-galaxy
        xx1<-numeric()
        for (i in 1:nrow(c1)){
          ng1<-as.matrix(which(gen[,1]==c1[i,1]))
          ng2<-gen[ng1,]
          ng3<-as.matrix(which(ng2[,2]==c1[i,2]))
          xx1<-rbind(xx1,ng2[ng3,])
        }
        xx1<-as.matrix(xx1)
        x1<-matrix(xx1[,3:(ncol(xx1))],,(ncol(xx1)-2))
        x1<-t(x1)
        x1<-as.matrix(x1)
        wow<-cbind(fx,x1)
        bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
        ef<-as.matrix(bbbb[(ncol(fx)+1):nrow(bbbb),1])
        y<-y-x1%*%ef
        
        if (flag==1)
        {code<-random(fx=fx,gen=gen,phe=y,kk=kk,CLO)
        }
        if (flag==0)
        {code<-fix(x=fx,gen=gen,y=y,kk=kk)
        }
        
        x0<-t(gen[,3:(ncol(gen))])
        bb<-code
        bb<-as.matrix(bb)
        aa<-numeric()
        for (i in 1:nrow(as.matrix(unique(gen[,1])))){
          
          mc<-which(gen[,1]==i)
          mc<-as.matrix(mc)
          for (j in 1:(nrow(mc)-2))
          {if (bb[mc[j+1],10]<bb[mc[j],10] & bb[mc[j+1],10]<bb[mc[j+2],10]){aa<-rbind(aa,mc[j+1])}
          }
          if (bb[mc[1],10]<bb[mc[2],10]){aa<-rbind(aa,mc[1])}
          if (bb[mc[nrow(mc)],10]<bb[mc[nrow(mc)-1],10]){aa<-rbind(aa,mc[nrow(mc)])}
        }
        
        mi<-code[aa,2:3]
        style<-numeric()
        for (i in 1:nrow(mi))
        {
          for (j in 1:nrow(xx1))
          {if (mi[i,1]==xx1[j,1] & mi[i,2]==xx1[j,2])
          {style<-rbind(style,aa[i])
          }
          }
        }
        aa<-as.matrix(setdiff(aa,style))
        xx<-x0[,aa]
        par<-ebayes_EM(fx,xx,y)
        cc<-as.matrix(which(par[,1]<=0.01))
        selectpos1<-which(par[,1]<=0.01)
        
        if (nrow(cc)>0)
        {
          name<-as.matrix(aa[cc,1])
          xxx<-as.matrix(x0[,name])
          y<-as.matrix(y)
          lod<-likelihood(fx,xxx,y)
          dd<-as.matrix(which(lod[,1]>=sLOD))
          if (nrow(dd)>0)
          {
            na<-as.matrix(name[dd])
            wow<-cbind(fx,xxx[,dd])
            bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
            ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
            galaxy2<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
            galaxyy<-rbind(galaxy,galaxy2)
            woww<-wow
          }else if (nrow(dd)==0)
          {galaxyy<-galaxy
          woww<-wow
          }
        }
        if (nrow(cc)==0)
        {galaxyy<-galaxy
        woww<-wow
        }
        pp<-as.matrix(phe[,1])
        va<-galaxyy[,3]*galaxyy[,3]
        ve<-(1/(n_sam-1))*t(pp-woww%*%bbbb)%*%(pp-woww%*%bbbb)
        vp<-(1/(n_sam-1))*t(pp-mean(pp))%*%(pp-mean(pp))
        vy<-(sum(va)+ve)
        
        if (vy>=vp){
          heredity<-va/as.vector(vy)
          pv<-vy}
        if (vy<vp){
          heredity<-va/as.vector(vp)
          pv<-vp}
        
        va<-matrix(va,,1)
        va<-round(va,4)
        heredity<-100*heredity
        heredity<-matrix(heredity,,1)
        heredity<-round(heredity,4)
        galaxyy[which(abs(galaxyy)>1e-4)]<-round(galaxyy[which(abs(galaxyy)>1e-4)],4)
        galaxyy[which(abs(galaxyy)<1e-4)]<-as.numeric(sprintf("%.4e",galaxyy[which(abs(galaxyy)<1e-4)]))
        
        if(is.null(mapname)==FALSE){
          map<-as.numeric(mapname[,2:3])
          map<-matrix(map,nrow(mapname),2)
          ###
          galaxytwo<-matrix(galaxyy[,1:2],,2)
          left_marker<-numeric()
          right_marker<-numeric()
          
          
          for( i in 1:nrow(galaxyy)){
            
            allchr<-as.vector(map[which(map[,1]==galaxytwo[i,1]),2])
            chr_loc<-which(map[,1]==galaxytwo[i,1])
            ###
            allmarker<-as.matrix(mapname[chr_loc,1])
            chose_left<-(map[which(map[,1]==galaxytwo[i,1]),2]<=galaxytwo[i,2])
            max_left<-max(allchr[chose_left[]==TRUE])
            chr_loclen<-length(chr_loc)
            if(max_left==-Inf){
              
              leftmarker<-matrix(mapname[chr_loc[1],1],,1)
            }else{
              leftloc<-which(allchr[]==max_left)
              
              leftmarker<-matrix(allmarker[leftloc],,1)
            }
            
            chose_right<-(map[which(map[,1]==galaxytwo[i,1]),2]>=galaxytwo[i,2])
            min_right<-min(allchr[chose_right[]==TRUE])
            if(min_right==Inf){
              rightmarker<-matrix(mapname[chr_loc[chr_loclen],1],,1)
            }else{
              rightloc<-which(allchr[]==min_right)
              
              rightmarker<-matrix(allmarker[rightloc],,1)
            }
            
            
            left_marker<-rbind(left_marker,leftmarker)
            right_marker<-rbind(right_marker,rightmarker)
          }
          
        }else{
          left_marker<-matrix("------",nrow(galaxyy),1)
          right_marker<-matrix("------",nrow(galaxyy),1) 
        }
        
        if((is.null(chrRaw_name)==FALSE)&&(is.null(chr_name)==FALSE)){
          chr_name<-chr_name
          chrRaw_name<-chrRaw_name
          galaxyysec<-galaxyy[,1]
          galaxyylast<-matrix(galaxyy[,2:ncol(galaxyy)],,(ncol(galaxyy)-1))
          chrName<-numeric()
          for( i in 1:length(galaxyysec)){
            chrLoc<-which(chr_name[]==galaxyysec[i])
            chrName0<-matrix(chrRaw_name[chrLoc],,1)
            chrName<-rbind(chrName,chrName0)
          }
          
          galaxyy_A<-cbind(chrName,galaxyylast)
        }else{
          galaxyy_A<-galaxyy
        }
        
        galaxyy_A<-as.matrix(galaxyy_A)
        vee<-matrix("",nrow(galaxyy_A),1)
        vee[1,1]<-round(ve,4)
        vee<-matrix(vee,,1)
        vpp<-matrix("",nrow(galaxyy_A),1)
        vpp[1,1]<-round(pv,4)
        vpp<-matrix(vpp,,1)
        traitid<-matrix(pheRaw[1,NUM+1],nrow(galaxyy_A),1)
        galaxyy<-galaxyy
        
        if(is.null(galaxyy)==FALSE){
          if(nrow(galaxyy)>1){
            galaxyy1<-galaxyy[,c(1,2,4)]
          }else{
            galaxyy1<-t(as.matrix(galaxyy[,c(1,2,4)]))
          }
        }
        
        result<-cbind(traitid,galaxyy_A,left_marker,right_marker,va,heredity,vee,vpp)
        colnames(result)<-c("Trait","Chr","Position (cM)","Additive Effect","LOD","Left_Marker","Right_Marker","Var_Genet_(i)","r2 (%)","Var_Error",
                            "Var_Phen (total)")
        
      }
    } 
    output<-list(result=result,mxmp=mxmp,galaxyy1=galaxyy1,res11=res11,chr_name=chr_name)
    return(output)
  }
}



