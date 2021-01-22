QTL.gCIMapping<-function(file=NULL,fileFormat="GCIM",fileICIMcov=NULL,Population=NULL,Model="Random",WalkSpeed=NULL,
               CriLOD=NULL,Likelihood="REML",SetSeed=11001,flagrqtl=FALSE,DrawPlot=TRUE,PlotFormat="jpeg",Resolution="Low",Trait=NULL,dir=NULL){
  
 
  WEN1re<-NULL; W1re<-NULL;readraw<-NULL;DoResult<-NULL;CLO<-NULL;
  
  readraw<-Readdata(file,fileFormat,fileICIMcov)
  DoResult<-Dodata(fileFormat,Population,Model,readraw)
  
  print("Running in progress, please be patient...")
  
  if(is.character(file)==FALSE){CLO<-1}
    
    pheRaw<-DoResult$pheRaw;genRaw<-DoResult$genRaw;mapRaw1<-DoResult$mapRaw1
    flag<-DoResult$flag;flagRIL<-DoResult$flagRIL;yygg1<-DoResult$yygg1;cov_en<-DoResult$cov_en
    
    if(Resolution=="Low"){
      widqqvalue<-1500
      heightqqvalue<-600
      pointsizeqqvalue<-12
      resppi<-72 
      ####################
      if(Population=="F2"){
        legend_size=1.2
        mainline_size=2.5
        backline_size=1.5
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }else{
        legend_size=1.2
        mainline_size=2.5
        backline_size=1.5
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }
    }else if(Resolution=="High"){
      widqqvalue<-10000
      heightqqvalue<-4000
      pointsizeqqvalue<-30
      resppi<-300 
      ########################### 
      if(Population=="F2"){
        legend_size=0.8
        mainline_size=2.5
        backline_size=0.8
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }else{
      legend_size=0.8
      mainline_size=2.5
      backline_size=2.5
      margin_space=1.5
      axis_space=1.0
      logPCoff=1.5
      lodthred=2.5
     }
    }
    gcimFunc <- function(mxmp,galaxyy1,res11,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,lodthred)
    {
      chr_pos <- mxmp[,1:2]
      chr_num <- length(chr_name)
      chr <- matrix(0,chr_num,1)
      pos <- matrix(0,chr_num,1)
      for(i in 1:chr_num)
      {
        temp <- numeric()
        temp <- length(which(chr_pos[,1]==i))
        if(i==1)
        {
          pos[i] <- temp
          chr[i] <- chr_pos[pos[i],2]
        }else{
          pos[i] <- pos[i-1] + temp
          chr[i] <- chr_pos[pos[i],2]
        }
      }
      ################# "pos" last ID, "pos" last value
      
      pos_acc <- matrix(0,chr_num,1)
      for(i in 1:chr_num)
      {
        if(i==1){
          pos_acc[i] <- chr[i]
        }else{
          pos_acc[i] <- pos_acc[i-1] + chr[i]
        }
      }
      
      firFil <- res11[,1:2]
      newposadd <- as.matrix(firFil[,2])
      for(i in 1:chr_num)
      {
        temp1 <- numeric()
        temp1 <- which(firFil[,1]==i)
        if(i>1)
        {
          newposadd[temp1] <- newposadd[temp1]+pos_acc[i-1]
        }
      }
      if(is.null(galaxyy1)==FALSE){
        if(is.null(dim(galaxyy1))==TRUE){
          galaxyy1<-matrix(galaxyy1,1,3)
        }
        newres_pos <- galaxyy1[,2]
        res_sumpos <- pos_acc[galaxyy1[which(galaxyy1[,1]>1),1]-1] + galaxyy1[which(galaxyy1[,1]>1),2]
        newres_pos[which(galaxyy1[,1]>1)] <- res_sumpos
        pospic<-c(newres_pos)
        lodpic<-c(galaxyy1[,3])
        resdf <- data.frame(pospic,lodpic)
      }
      
      resp<-as.matrix(res11[,3])
      pmin<-min(resp[resp!=0])
      locsub<-which(resp==0)
      if(length(locsub)!=0){
        subvalue<-10^(1.1*log10(pmin))
        res11[locsub,3]<-subvalue
      }else{
        res11<-res11
      }
      
      negloP <- -log10(as.matrix(res11[,3]))
      
      if(is.null(galaxyy1)==FALSE){
        ###################################change 20200914 two y axis
        par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
        plot(newposadd,negloP,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
             xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
        axis(side=4,cex.axis=legend_size)
        abline(v=pos_acc,lty=2,col="gray")
        par(new=TRUE)
        plot(pospic,lodpic,type="h",col=color1,yaxt="n",xlab="Genome position (cM)",ylab="",
             cex.axis=legend_size,cex.lab=legend_size+0.4,lwd=mainline_size,xlim=c(0,max(newposadd)),
             ylim=c(0,max(lodpic)))
        axis(side=2,cex.axis=legend_size)
        mtext("LOD score",side=2,line=3*axis_space,cex=legend_size+0.4,col=color1)
        abline(h=lodthred,lty=5,col=color2)
        mtext(expression(-'log'[10]*'(P-value)'),side=4,line=3*axis_space,cex=legend_size+0.4,col = color2)
        
      }else{
        par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
        plot(newposadd,negloP,type="l",col=color1,xlab="Genome position (cM)",ylab=expression('Expected -log'[10]*'(P)'),cex.axis=legend_size+0.4,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
        abline(h=lodthred,lty=5,col=color1)
      }
    }
    
    
    gcimFuncF2 <- function(mxmp,galaxyy1,res1a,res1d,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,color3,lodthred)
    {
    chr_pos <- mxmp[,1:2]
    chr_num <- length(chr_name)
    chr <- matrix(0,chr_num,1)
    pos <- matrix(0,chr_num,1)
    for(i in 1:chr_num)
    {
      temp <- numeric()
      temp <- length(which(chr_pos[,1]==i))
      if(i==1)
      {
        pos[i] <- temp
        chr[i] <- chr_pos[pos[i],2]
      }else{
        pos[i] <- pos[i-1] + temp
        chr[i] <- chr_pos[pos[i],2]
      }
    }
    
    pos_acc <- matrix(0,chr_num,1)
    for(i in 1:chr_num)
    {
      if(i==1){
        pos_acc[i] <- chr[i]
      }else{
        pos_acc[i] <- pos_acc[i-1] + chr[i]
      }
    }
    firFila <- res1a[,1:2]
    newposadda <- as.matrix(firFila[,2])
    for(i in 1:chr_num)
    {
      temp1a <- numeric()
      temp1a <- which(firFila[,1]==i)
      if(i>1)
      {
        newposadda[temp1a] <- newposadda[temp1a]+pos_acc[i-1]
      }
    }
    # firFild <- res1d[,1:2]
    # newposaddd <- as.matrix(firFild[,2])
    # for(i in 1:chr_num)
    # {
    #   temp1d <- numeric()
    #   temp1d <- which(firFild[,1]==i)
    #   if(i>1)
    #   {
    #     newposaddd[temp1d] <- newposaddd[temp1d]+pos_acc[i-1]
    #   }
    # }
    newposaddd<-newposadda
    #############newposaddd==newposadda
    if(is.null(galaxyy1)==FALSE){
      if(is.null(dim(galaxyy1))==TRUE){
        galaxyy1<-matrix(galaxyy1,1,3)
      }
      newres_pos <- galaxyy1[,2]
      res_sumpos <- pos_acc[galaxyy1[which(galaxyy1[,1]>1),1]-1] + galaxyy1[which(galaxyy1[,1]>1),2]
      newres_pos[which(galaxyy1[,1]>1)] <- res_sumpos
      pospic<-c(newres_pos)
      lodpic<-c(galaxyy1[,3])
      resdf <- data.frame(pospic,lodpic)
    }
    negloPa <- as.matrix(res1a[,3])
    negloPd <- as.matrix(res1d[,3])
    if(is.null(galaxyy1)==FALSE){
      #####################################change 20200914 two y axis
      par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
      plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
           xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
      par(new=TRUE)
      plot(newposaddd,negloPd,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
           xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
      axis(side=4,at=seq(0,logPCoff*max(negloPa,negloPd),ceiling(logPCoff*max(negloPa,negloPd)/5)),cex.axis=legend_size)
      abline(v=pos_acc,lty=2,col="gray")
      par(new=TRUE)
      plot(pospic,lodpic,type="h",col=color1,yaxt="n",xlab="Genome position (cM)",ylab="",
           cex.axis=legend_size,cex.lab=legend_size+0.4,lwd=mainline_size,xlim=c(0,max(newposadda,newposaddd)),
           ylim=c(0,max(lodpic)))
      axis(side=2,cex.axis=legend_size)
      mtext("LOD score",side=2,line=3*axis_space,cex=legend_size+0.4,col=color1)
      abline(h=lodthred,lty=5,col=color2)
      mtext(expression(-'log'[10]*'(P-value)'),side=4,line=3*axis_space,cex=legend_size+0.4,col = color2)
    }else{
      ##########################change 20200914
      par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
      plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
           xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
      par(new=TRUE)
      plot(newposaddd,negloPd,type="l",col=color2,yaxt="n",xlab="Genome position (cM)",ylab="",lwd=backline_size,
           xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
      axis(side=4,at=seq(0,logPCoff*max(negloPa,negloPd),ceiling(logPCoff*max(negloPa,negloPd)/5)),cex.axis=legend_size)
      abline(v=pos_acc,lty=2,col="gray")
      mtext(expression(-'log'[10]*'(P-value)'),side=4,line=3*axis_space,cex=legend_size+0.4,col = color2)
      
    }
  }
    
    
    if(Population=="F2"){
      WEN1re<-WenF(pheRaw,genRaw,mapRaw1,yygg1,cov_en,WalkSpeed,CriLOD,dir)
      for(NUM in Trait){
        rewen<-NULL;mxmp=NULL;galaxyy1<-NULL;res1a=NULL;res1d=NULL;chr_name=NULL 
        TRY1<-try({
          outWEN<-WenS(flag,CriLOD,NUM,pheRaw,Likelihood,SetSeed,flagrqtl,WEN1re$yygg,WEN1re$mx,WEN1re$phe,WEN1re$chr_name,
                       WEN1re$v.map,WEN1re$gen.raw,WEN1re$a.gen.orig,WEN1re$d.gen.orig,WEN1re$n,WEN1re$names.insert2,WEN1re$X.ad.tran.data,WEN1re$X.ad.t4,dir)
          rewen<-outWEN$result
          mxmp<-outWEN$mxmp;galaxyy1<-outWEN$galaxyy1;res1a<-outWEN$res1a;res1d<-outWEN$res1d;chr_name<-outWEN$chr_name
          
        },silent=FALSE)    
        
        if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){  
          TRY2<-try({ 
              
            write.table(rewen,paste(dir,"/",NUM,"_GCIM result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
          
            if(DrawPlot==TRUE){
              if(PlotFormat=="png")
              {
                png(paste(dir,"/",NUM,"_resF2.png",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="tiff"){
                tiff(paste(dir,"/",NUM,"_resF2.tiff",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="jpeg"){
                jpeg(paste(dir,"/",NUM,"_resF2.jpeg",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="pdf"){
                pdf(paste(dir,"/",NUM,"_resF2.pdf",sep=""), width=16)
              }
              gcimFuncF2(mxmp,galaxyy1,res1a,res1d,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,"red","gray50","green",lodthred)
               
              dev.off()
            }
          },silent=FALSE)  
        }
      }
      
    }else{
      
      W1re<-WangF(pheRaw,genRaw,mapRaw1,yygg1,flagRIL,cov_en,Population,WalkSpeed,CriLOD)
      
      for(NUM in Trait){
        rew<-NULL;mxmp=NULL;galaxyy1<-NULL;res11=NULL;chr_name=NULL
        TRY1<-try({
          outW<-WangS(flag,CriLOD,NUM,pheRaw,W1re$chrRaw_name,W1re$yygg,W1re$mx,W1re$phe,W1re$chr_name,W1re$gen,W1re$mapname,CLO)
          rew<-outW$result
          mxmp<-outW$mxmp;galaxyy1<-outW$galaxyy1;res11<-outW$res11;chr_name<-outW$chr_name
        },silent=FALSE)  
        
        
        if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){   
          TRY2<-try({ 
            write.table(rew,paste(dir,"/",NUM,"_GCIM result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
            
            if(DrawPlot==TRUE){
              if(PlotFormat=="png")
              {
                png(paste(dir,"/",NUM,"_res.png",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="tiff"){
                tiff(paste(dir,"/",NUM,"_res.tiff",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="jpeg"){
                jpeg(paste(dir,"/",NUM,"_res.jpeg",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="pdf"){
                pdf(paste(dir,"/",NUM,"_res.pdf",sep=""), width=16)
              }
              
              gcimFunc(mxmp,galaxyy1,res11,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,"red","gray50",lodthred)
 
              dev.off()
            }
          },silent=FALSE)  
       }
    }  
  }
}
