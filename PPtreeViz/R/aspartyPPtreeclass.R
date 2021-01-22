
#' @export
#' @import partykit 
#' 
as.party.PPtreeclass<-function(obj,data=TRUE,Rule=1,...){
  ff <- data.frame(obj$Tree.Struct)
  n <- nrow(ff)
  if (n == 1) 
    return(partykit::partynode(as.integer(1)))
  is.leaf <- (ff$Index == 0)
  ncomplete<-rep(2,n)
  ncomplete[is.leaf]<-0
  index<-cumsum(c(1,ncomplete+1*(!is.leaf)))
  primary<-numeric(n)
  primary[!is.leaf]<-index[c(!is.leaf,FALSE)]
  mf<-obj$origdata%*%t(obj$projbest.node)
  rownames(mf)<-rownames(obj$origdata)
  colnames(mf)<-paste("proj",1:ncol(mf),sep="")
  mf<-data.frame(mf)
  PPtreeclass_fitted <- function() {
    ret <- as.data.frame(matrix(nrow = NROW(mf), ncol = 0))
    ret[["(fitted)"]] <- apply(matrix(as.numeric(predict.PPtreeclass(obj)),ncol=1),1,
                               function(x) which((ff$R.F.node.ID==x)*is.leaf==1))
    ret[["(response)"]] <- obj$origclass
    ret
  }
  
  fitted <- PPtreeclass_fitted()
  PPtreeclass_kids <- function(i) {
    if (is.leaf[i]) 
      return(NULL)
    else 
      return(c(ff[i,c(3,2)])) 
  }
  
  PPtreeclass_split <- function(j) {
    if (j < 1) 
      return(NULL)
    idj <- as.integer(ff$Coef.ID[j])
    ret <- partykit::partysplit(varid = idj, 
                      breaks = as.double(obj$splitCutoff.node[idj,Rule]), 
                      right = FALSE, 
                      index = 2L:1L)
    ret
  }
  
  PPtreeclass_node <- function(i) {
    if (is.null(PPtreeclass_kids(i))) 
      return(partynode(as.integer(i)))
    nd <- partykit::partynode(as.integer(i), split = PPtreeclass_split(i), 
                    kids = lapply(PPtreeclass_kids(i),PPtreeclass_node))
    
    left <- partykit::nodeids(kids_node(nd)[[1L]], terminal = TRUE)
    right <- partykit::nodeids(kids_node(nd)[[2L]], terminal = TRUE)
    nd$split$prob <- c(0, 0)
    nl <- sum(fitted[["(fitted)"]] %in% left)
    nr <- sum(fitted[["(fitted)"]] %in% right)
    if(nl > nr) {
      nd$split$prob <- c(1, 0)
    } else {
      nd$split$prob <- c(0, 1)
    }
    nd$split$prob <- as.double(nd$split$prob)
    return(nd)
  }
  node <- PPtreeclass_node(1)
  rval <- partykit::party(node = node, 
                data = if (data) 
                  mf
                else mf[0L, ], 
                fitted = fitted, 
                terms = obj$terms, 
                info = list(method = "PPtreeclass"))
  class(rval) <- c("constparty", class(rval))
  rval
  return(rval)
}