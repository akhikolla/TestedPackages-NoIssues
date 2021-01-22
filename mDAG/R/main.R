MatrixtoGraph_undirected=function(skeleton,data){
  graph=list()
  arcs <- which(skeleton!=0,arr.ind = T)
  arcs=t(apply(arcs,1,function(z) colnames(data)[z]))
  colnames(arcs)=c('from','to')
  rownames(arcs)=NULL
  graph$arcs=arcs
  for (i in 1:ncol(data)){
    nbr=arcs[which(arcs[,2]==colnames(data)[i]),1]
    graph$nodes[[colnames(data)[i]]]$nbr=nbr
  }
  colnames(skeleton)=rownames(skeleton)=colnames(data)
  graph$skeleton=skeleton
  return(graph)
}

MatrixtoGraph_directed=function(skeleton,data){
  graph=list()
  arcs <- which(skeleton!=0,arr.ind = T)
  arcs=t(apply(arcs,1,function(z) colnames(data)[z]))
  colnames(arcs)=c('from','to')
  rownames(arcs)=NULL
  graph$arcs=arcs
  for (i in 1:ncol(data)){
    parents=arcs[which(arcs[,2]==colnames(data)[i]),1]
    children=arcs[which(arcs[,1]==colnames(data)[i]),2]
    
    graph$nodes[[colnames(data)[i]]]$parents=parents
    graph$nodes[[colnames(data)[i]]]$children=children
    graph$nodes[[colnames(data)[i]]]$nbr=union(parents,children)
  }
  colnames(skeleton)=rownames(skeleton)=colnames(data)
  graph$skeleton=skeleton
  return(graph)
}


mgm_skeleton=function(data = data, type = type, level = level, SNP = SNP ,lambdaGam = lambdaGam, 
                      lambdaSel = lambdaSel, ruleReg=ruleReg, alphaSel = alphaSel, threshold=threshold, weights){
  
  edgeWeights=pair_fast_mgm(data = data,
                            type = type,
                            level = level,lambdaGam = lambdaGam,
                            lambdaSel = lambdaSel,ruleReg=ruleReg,alphaSel = alphaSel,threshold=threshold )
  
  edgeWeights=(edgeWeights>0)
  edgeWeights[edgeWeights==T]=1
  edgeWeights[edgeWeights==F]=0
  
  # we don't allow edges among SNPs 
  SNP_ind=which(SNP==1)
  edgeWeights[SNP_ind,SNP_ind]=0
  return (MatrixtoGraph_undirected(edgeWeights,data))
  
}

penpc_skeleton=function(data,type,level,edgeWeights,indepTest=ConditionalTestPermute,nperm=nperm,alpha=alpha){
  
  suffStat=list('dat_in'=data,'type'=type,'level'=level,C = cor(data), n = nrow(data),'nperm'=nperm)
  
  penpcskel=skeletonPENstable_modified(suffStat, indepTest, p=as.integer(ncol(data)), alpha=alpha,
                                       edgeWeights, verbose=F)
  
  skeleton_penpc=wgtMatrix(penpcskel@graph)
  
  return (list("graph"=MatrixtoGraph_undirected(skeleton_penpc,data),"penpcskel"=penpcskel))
  
}


greedysearch_orientation=function(data,type,level,SNP,result,weights=rep(1, nrow(data))){
  
  estimateundirect=which(result!=0,arr.ind = T)
  arcs=t(apply(estimateundirect,1,function(z) colnames(data)[z]))
  colnames(arcs)=c('from','to')
  rownames(arcs)=NULL
  find_nbr=sapply(1:ncol(data),function(i) 
    colnames(data)[(estimateundirect[which(estimateundirect[,1]==i),2])] )
  
  
  rst=list('learning'=NULL,'nodes'=list(),'arcs'=arcs)
  rst$learning=list('whitelist'=NULL,'blacklist'=NULL,'test'='cor','optimized'=F,'ntests'=NULL,'algo'='mgmPenPC')
  
  for(i in 1:ncol(data)){
    mb=find_nbr[[i]]
    nbr=find_nbr[[i]]
    nbr_index=(1:ncol(data))[(colnames(data)%in%nbr)]
    parents=character(0)
    children=character(0)
    rst$nodes[[colnames(data)[i]]]=list('mb'=mb,'nbr'=nbr,'nbr_index'=nbr_index,'parents'=parents,'children'=children)
  }
  
  rst$skeleton=result
  
  
  skeleton=GreedySearch(data ,type, level, SNP, rst, weights)
  return(MatrixtoGraph_directed(skeleton,data))
}
