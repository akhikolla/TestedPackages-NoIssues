Readdata<-function(file=NULL,fileFormat=NULL,fileICIMcov=NULL){
  geoo=NULL;pho=NULL;poss=NULL;parm=NULL;y_jun3=NULL;genRaw1=NULL;pheRaw=NULL;mapRaw11=NULL;cov_en=NULL
  if(fileFormat=="ICIM"){
    geoo<-as.matrix(read.xlsx(file,sheet = "Genotype",colNames = FALSE))
    pho<-as.matrix(read.xlsx(file,sheet = "Phenotype",colNames = FALSE))
    poss<-read.xlsx(file,sheet = "LinkageMap",colNames = FALSE)
    parm<-read.xlsx(file,sheet = "GeneralInfo",colNames = FALSE)
    if(is.null(fileICIMcov)==FALSE){
      cov_en<-fread(fileICIMcov,header = FALSE,stringsAsFactors=T)
      cov_en<-as.matrix(cov_en)
    }
  }else if(fileFormat=="Cart"){
    y_jun3<-scan(file,what = "",sep = "")
    
  }else if(fileFormat=="GCIM"){
    if(is.character(file)==TRUE){
      genRaw<-fread(file,header = FALSE,stringsAsFactors=T)
    }else{
      genRaw<-file
    }
    
    titlenameGen<-genRaw[1,1:3]
    hapName<-c("marker","chr","pos")
    
    if(all(titlenameGen==hapName)==FALSE){
      warning("please check the Linkage map's name in file!")
    }
    traitloc<-which(genRaw[,2]=="trait1")
    
    if(length(traitloc)==0){
      warning("please check the phenotype in file!") 
    }
    envirloc<-which(genRaw[,2]=="Covar1")
    if(length(envirloc)!=0){
      pheRaw<-t(rbind(genRaw[1,][,-(1:2)],genRaw[traitloc:(envirloc-1),][,-(1:2)]))
      cov_en<-cbind(t(genRaw[1,][,-(1:2)]),t(genRaw[envirloc:nrow(genRaw),][,-(1:2)]))
      colnames(cov_en)<-NULL;rownames(cov_en)<-NULL
    }else{
      pheRaw<-t(rbind(genRaw[1,][,-(1:2)],genRaw[traitloc:nrow(genRaw),][,-(1:2)]))  
    }
    colnames(pheRaw)<-NULL;rownames(pheRaw)<-NULL
    genRaw1<-as.matrix(genRaw[1:(traitloc-1),-c(2,3)])
    mapRaw11<-as.matrix(genRaw[1:(traitloc-1),1:3])
  }
  result<-list(geoo=geoo,pho=pho,poss=poss,parm=parm,y_jun3=y_jun3,genRaw1=genRaw1,pheRaw=pheRaw,mapRaw11=mapRaw11,cov_en=cov_en)  
  return(result)
}

