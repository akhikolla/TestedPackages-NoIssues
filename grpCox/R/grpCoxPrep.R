#####################################################
#Information needed for Cox model and group penalty##
#####################################################
sglcoxinfo<-function(x, y, g, m, standardize = TRUE){
  
  # Extract info for Cox model
  ox <- x
  N0 <- nrow(x)
  # oi <- order(y[, "status"], decreasing=TRUE)
  # x <- x[oi, ]
  # y <- y[oi, ]
  oi <- order(y[, "time"])
  x <- x[oi, ]
  y <- y[oi, ]
  
  # Remove the first censored cases
  i1 <- which(y[, "status"]==1)
  if(length(i1)==0){
    x <- x; y <- y
  } else {
    mi1 <- min(i1)-1
    if (mi1!=0) {
      x <- x[-c(1:mi1), ];y <- y[-c(1:mi1), ]
    }
  }
  ty <- y[, "time"];tevent <- y[, "status"]
  N <- nrow(x);n1 <- sum(y[, "status"])
  
  dty <- duplicated(ty) # ties
  
  # For calculation of log-likelihood
  if (any(dty)) {
    tevent0 <- tevent
    tevent0[which(dty)] <- 0
    ievent <- cumsum(tevent0);loc1 <- which(tevent0==1)
    nevent <- table(ievent);n <- length(unique(ievent))
    nevent1 <- tapply(tevent==1, ievent, sum)
  } else {
    ievent <- cumsum(tevent);loc1 <- which(tevent==1)
    nevent <- table(ievent);n <- length(unique(ievent))
    nevent1 <- rep(1, n)
  }
  x  <-  as.matrix(x)
  
  # Feature-level standardization
  if(standardize == TRUE){
    tem  <-  stdQ(x)
    XX <- tem$x
    scale <- tem$sd
    
    # Extract info for group penalty
    # Setup group
    G <- setupG(g, m)
    
    # Re-order groups if necessary
    G  <-  reorderG(G, attr(G, 'm'))
    if(attr(G,'reorder')) XX  <-  XX[,attr(G,'ord')]
    
    # Set group multiplier if missing
    m <- attr(G, 'm')
    if (all(is.na(m))) {
      m <- sqrt(table(g[g!=0])) 
    }
    
    K  <-  as.numeric(table(g))
    K1  <-  as.integer(c(0,cumsum(K)))
    p <- dim(XX)[2]
    
    return(list(x=XX, N0=N0, tevent=tevent, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n, 
                g=g, m=m, K=K, K1=K1, ox=ox, p=p,
                reorder=attr(G,'reorder'), ord.inv=attr(G,'ord.inv'), scale=scale))
  }else{
    XX <- x
    # Extract info for group penalty
    # Setup group
    G <- setupG(g, m)
    
    # Re-order groups if necessary
    G  <-  reorderG(G, attr(G, 'm'))
    if(attr(G,'reorder')) XX  <-  XX[,attr(G,'ord')]
    
    # Set group multiplier if missing
    m <- attr(G, 'm')
    if (all(is.na(m))) {
      m <- sqrt(table(g[g!=0])) 
    }
    
    K  <-  as.numeric(table(g))
    K1  <-  as.integer(c(0,cumsum(K)))
    p <- dim(XX)[2]
    return(list(x=XX, N0=N0, tevent=tevent, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n, 
                g=g, m=m, K=K, K1=K1, ox=ox, p=p,  
                reorder=attr(G,'reorder'), ord.inv=attr(G,'ord.inv')))
  }
}

setupG <- function(group, m) {
  gf <- factor(group)
  g <- as.integer(gf)
  lev <- levels(gf)
  
 if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups")
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
  }
  structure(g, levels=lev, m=m)
}

reorderG <- function(g, m) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    g <- as.integer(gf)
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels=lev, m=m, ord=ord, ord.inv=ord.inv, reorder=reorder)
}


#####################################################
#  Overlapping processing procedure                ##
#####################################################
### expanded latent matrix 
expandedlatent <- function(x, group){
  # check whether or not group is a list
  if(class(group) != 'list'){
    stop("group must be a list of numeric indices of variables")
  }
  # number of group
  J <- length(group)
  # number of variables
  P <- ncol(x)
  
  # variables and groups information 
  var_grp <- Matrix(0, nrow = J, ncol = P, sparse = TRUE, 
                   dimnames = list(as.character(rep(NA,J)), as.character(rep(NA,P))))
  colnames(x) <- paste('X', 1:P, sep = "")
  if(is.null(names(group))){
    names(group) <- paste('grp', 1:J, sep = "")
  }
  for(j in 1:J){
    ind <- group[[j]]
    var_grp[j, ind] <- 1
    colnames(var_grp)[ind] <- colnames(x)[ind]
  }
  rownames(var_grp) <- as.character(names(group))
  
  # convert group into group-vector type
  # number of variables in each group
  num_vargrp <- apply(var_grp, 1, sum)
  grp_vec <- rep(1:J, times=num_vargrp)
  
  # expanded matrix x_expand
  x <- as.matrix(x)
  x_expand <- NULL
  names <- NULL
  for (i in 1:J){
    idx = var_grp[i,]==1
    x_expand <- cbind(x_expand, x[,idx,drop=FALSE])
    names <- c(names, colnames(var_grp)[idx])
  }
  colnames(x_expand) <- paste('grp', grp_vec, '_', names=names, sep = "")
  attr(x_expand, 'group') <- group
  attr(x_expand, 'var_grp') <- var_grp 
  attr(x_expand, 'grp_vec') <- grp_vec
  structure(x_expand, group=group, var_grp=var_grp, grp_vec=grp_vec)
}

### convert the latent coefficients in expanded space back to original coefficients
latent2ori <- function(coef.latent, var_grp, grp_vec){
  P <- ncol(var_grp)
  J <- nrow(var_grp)
  
  coef.ori <- matrix(0, ncol = ncol(coef.latent), nrow = P)
  for(j in 1:J){
    ind <- which(var_grp[j,] == 1)
    coef.ori[ind,] <- coef.ori[ind,] + coef.latent[which(grp_vec == j),,drop=FALSE]
  }
  rownames(coef.ori) <- colnames(var_grp)
  coef.ori
}


#############################################
#####  Plot coefficients' paths         #####
#############################################
plot.Coef <- function(x,lambda,label=TRUE,xlab="log(Lambda)",ylab="Coefficients",title=NULL,...){
  which <- nonzeroCoef(x)
  nwhich<-length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={
           warning("No plot produced since all coefficients zero")
           return()
         },
         "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  x<-as.matrix(x[which,,drop=FALSE])
  index<-log(lambda)
  # library(colorspace)
  col_set <- rainbow_hcl(nrow(x))
  matplot(index,t(x),
          type = "l", lty = "solid", lwd = 3, ylab = "Coefficients",
          xlab = "log(Lambda)", main = title, col = col_set,
          cex.axis = 1.25, font = 2, cex.lab = 1.5, col.lab = '#993333', font.lab=2)
  
  if(label){
    nnz<-length(which)
    xpos<-min(index)
    pos<-3
    xpos<-rep(xpos,nnz)
    ypos<-x[,ncol(x)]
    text(xpos,ypos-.05,paste(which),cex=.8,pos=pos)
  }
}

plot.gCoef <- function(x,g,lambda,label=TRUE,xlab="log(Lambda)",ylab="Coefficients",title=NULL,...){
  which <- nonzeroCoef(x)
  nwhich<-length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={
           warning("No plot produced since all coefficients zero")
           return()
         },
         "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  x<-as.matrix(x[which,,drop=FALSE])
  index<-log(lambda)
  # library(colorspace)
  set <- rainbow_hcl(length(unique(g)))
  col_set <- rep(0, length(g))
  for(i in 1:length(g)){
    col_set[i] <- set[g[i]]
  }
  
  matplot(index,t(x),
          type = "l", lty = "solid", lwd = 2, ylab = "Coefficients",
          xlab = "log(Lambda)", main = title, col = col_set,
          cex.axis = 1.25, font = 2, cex.lab = 1.5, col.lab = '#993333', font.lab=2)
  
  if(label){
    nnz<-length(which)
    xpos<-min(index)
    pos<-3
    xpos<-rep(xpos,nnz)
    ypos<-x[,ncol(x)]
    text(xpos,ypos-.05,paste(which),cex=.8,pos=pos)
  }
}

nonzeroCoef <- function (beta, bystep = FALSE) 
{
  ### bystep = FALSE means which variables were ever nonzero
  ### bystep = TRUE means which variables are nonzero for each step
  nr<-nrow(beta)
  if (nr == 1) {#degenerate case
    if (bystep) 
      apply(beta, 2, function(x) if (abs(x) > 0) 
        1
        else NULL)
    else {
      if (any(abs(beta) > 0)) 
        1
      else NULL
    }
  }
  else {
    beta<-abs(beta)>0 # this is sparse
    which<-seq(nr)
    ones<-rep(1,ncol(beta))
    nz<-as.vector((beta%*%ones)>0)
    which<-which[nz]
    if (bystep) {
      if(length(which)>0){
        beta<-as.matrix(beta[which,,drop=FALSE])
        nzel <- function(x, which) if (any(x)) 
          which[x]
        else NULL
        which<-apply(beta, 2, nzel, which)
        if(!is.list(which))which<-data.frame(which)# apply can return a matrix!!
        which
      }
      else{
        dn<-dimnames(beta)[[2]]
        which<-vector("list",length(dn))
        names(which)<-dn
        which
      }
    }
    else which
  }
}

#############################################
#####  Plot log likelihood paths        #####
#############################################
plot.llCV = function(x,...){
  error.bars =function(object, upper, lower, width = 0.02, ...)
  {
    xlim = range(object)
    barw = diff(xlim) * width
    segments(object, upper, object, lower, ...)
    segments(object - barw, upper, object + barw, upper, ...)
    segments(object - barw, lower, object + barw, lower, ...)
    range(upper, lower)
  }
  
  res = x$fit
  plot.args=list(x=log(res[,1]),y=res[,2],ylim=range((res[,2]+res[,3]),(res[,2]-res[,3])),
                 xlab="log(Lambda)",ylab="CV Partial log-likelihood",type="n",
                 cex.axis = 1.25, font = 2, cex.lab = 1.5, col.lab = '#993333', font.lab=2)
  do.call("plot",plot.args)
  error.bars(log(res[,1]),res[,2]+res[,3],res[,2]-res[,3],width=0.01,col="darkgrey")
  points(log(res[,1]),res[,2],pch=20,col="red")
  lambda_pcvl = which(res[,4]=="pcvl")
  if(length(lambda_pcvl)==0){
    lambda_pcvl = which(res[,4]=="pcvl=cvmax")
    lambda_max = which(res[,4]=="pcvl=cvmax")
  }else{
    lambda_max = which(res[,4]=="cvmax")
  }
  abline(v=log(res[,1][lambda_max]),lty=3)
  abline(v=log(res[,1][lambda_pcvl]),lty=3)
  invisible()
}