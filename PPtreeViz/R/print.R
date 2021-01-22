#' Print PP.Tree.class result
#' 
#' Print the projection pursuit classification tree result
#' @title Print PP.Tree.class result
#' @param x PPtreeclass object
#' @param coef.print print projection coefficients in each node ifTRUE
#' @param cutoff.print print cutoff values in each node if TRUE
#' @param verbose print if TRUE, no output if FALSE
#' @param ... arguments to be passed to methods
#' @references Lee, EK(2017) 
#' PPtreeViz: An R Package for Visualizing Projection Pursuit Classification 
#' Trees, Journal of Statistical Software <doi:10.18637/jss.v083.i08>
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- PPTreeclass(Species~.,data = iris,"LDA")
#' Tree.result
#' print(Tree.result,coef.print=TRUE,cutoff.print=TRUE)
#' 
print.PPtreeclass<-function(x,coef.print=FALSE,cutoff.print=FALSE,
                            verbose=TRUE,...){
  PPtreeOBJ<-x
  TS<-PPtreeOBJ$Tree.Struct
  Alpha<-PPtreeOBJ$projbest.node
  cut.off<-PPtreeOBJ$splitCutoff.node
  gName<-names(table(PPtreeOBJ$origclass))
  pastemake<-function(k,arg,sep.arg=""){
    temp<-""
    for(i in 1:k)
      temp<-paste(temp,arg,sep=sep.arg)
    return(temp)
  }
  TreePrint<-"1) root"
  i<-1  
  flag.L<-rep(FALSE,nrow(TS))
  keep.track<-1
  depth.track<-0
  depth<-0
  while(sum(flag.L)!=nrow(TS)){
    if(!flag.L[i]){                    
      if(TS[i,2] == 0) {
        flag.L[i]<-TRUE
        n.temp<-length(TreePrint)
        tempp<-strsplit(TreePrint[n.temp],") ")[[1]]
        temp.L<-paste(tempp[1],")*",tempp[2],sep="")
        temp.L<- paste(temp.L,"  ->  ","\"",gName[TS[i,3]],"\"",sep="")
        TreePrint<-TreePrint[-n.temp]
        id.l<-length(keep.track)-1
        i<-keep.track[id.l]
        depth<-depth -1
      } else if(!flag.L[TS[i,2]]){
        depth<-depth+1
        emptyspace<-pastemake(depth,"   ")
        temp.L<-paste(emptyspace,TS[i,2],")  proj",
                      TS[i,4],"*X < cut",TS[i,4],sep="")
        i<-TS[TS[i,2],1]   
      } else{
        depth<-depth +1
        emptyspace<-pastemake(depth,"   ")          
        temp.L<- paste(emptyspace,TS[i,3],")  proj",
                       TS[i,4],"*X >= cut",TS[i,4],sep="")
        flag.L[i]<-TRUE
        i<-TS[TS[i,3],1]
      } 
      keep.track<-c(keep.track,i)
      depth.track<-c(depth.track,depth)
      TreePrint<-c(TreePrint,temp.L)
    } else{
      id.l<-id.l-1
      i<-keep.track[id.l]
      depth<-depth.track[id.l]
    }
  }
  colnames(Alpha)<-colnames(PPtreeOBJ$origdata)
  rownames(Alpha)<-paste("proj",1:nrow(Alpha),sep="")  
  colnames(cut.off)<-paste("Rule",1:ncol(cut.off),sep="")
  rownames(cut.off)<-paste("cut",1:nrow(cut.off),sep="")
  TreePrint.output<-
    paste("=============================================================",                          
          "\nProjection Pursuit Classification Tree result",                           
          "\n=============================================================\n")
  for(i in 1:length(TreePrint))
    TreePrint.output<-paste(TreePrint.output,TreePrint[i],sep="\n")
  TreePrint.output<-paste(TreePrint.output,"\n",sep="")
  sample.data.X<-PPtreeOBJ$origdata
  sample.data.class<-PPtreeOBJ$origclass
  error.rate<-matrix(c(PPclassify(PPtreeOBJ,sample.data.X, Rule=1, 
                                  sample.data.class)$predict.error,
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=2, 
                                  sample.data.class)$predict.error,
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=3, 
                                  sample.data.class)$predict.error,
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=4, 
                                  sample.data.class)$predict.error,
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=5, 
                                  sample.data.class)$predict.error,
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=6, 
                                  sample.data.class)$predict.error,
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=7, 
                                  sample.data.class)$predict.error,                 
                       PPclassify(PPtreeOBJ,sample.data.X, Rule=8, 
                                  sample.data.class)$predict.error)/
                       nrow(sample.data.X),nrow=1)               
  colnames(error.rate)<-colnames(cut.off)
  rownames(error.rate)<-c("error.rate")
  colnames(Alpha)<-paste(1:ncol(Alpha),":\"",colnames(Alpha),"\"",sep="")
  if(verbose){
    cat(TreePrint.output)
    if(coef.print){
      cat("\nProjection Coefficient in each node",
          "\n-------------------------------------------------------------\n")
      print(round(Alpha,4))
    }
    if(cutoff.print){
      cat("\nCutoff values of each node",
          "\n-------------------------------------------------------------\n")
      print(round(cut.off,4))
    }    
    cat("\nError rates of various cutoff values",
        "\n-------------------------------------------------------------\n")
    print(round(error.rate,4))
  }    
  return(invisible(TreePrint)) 
}
#' @export
print.PPoptim<-function(x,index.print=TRUE,proj.print=TRUE,orig.print=FALSE,...){
  PPoptOBJ<-list()
  if(index.print){
   PPoptOBJ<-c(PPoptOBJ,list(indexbest=x[[1]]))
  } 
  if(proj.print){
   PPoptOBJ<-c(PPoptOBJ,list(projbest=x[[2]]))
  } 
  if(orig.print){
   PPoptOBJ<-c(PPoptOBJ,list(origclass=x[[3]],origdata=x[[4]]))
  } 
  print(PPoptOBJ)
}

   