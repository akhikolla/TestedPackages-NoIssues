## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
PACE2d = function(mat, ev.vec, sigma2, ev.val, var.est=FALSE){
  i = NULL
  res = foreach (i = 1:nrow(mat)) %dopar%{
    nonna.idx = !is.na(mat[i,])
    non.na.num = sum(nonna.idx)
    ## could be a problem when sigma2 is small
    ## sigma2 = max(0.5, sigma2)
    sigma2 = max(0.0001, sigma2)
    phi = ev.vec[nonna.idx,, drop=FALSE]
    H = ev.val * t(phi)
    sigma = phi %*% H + diag(sigma2, non.na.num)
    list(
      ## estimation of xi
      H %*% solve(sigma, mat[i,nonna.idx]),
      ## variance estimation of xi
      H %*% solve(sigma, t(H))
    )
  }
  ##estimated xi matrix: K x length(idseval)
  xi.est = do.call('cbind', lapply(1:length(res), function(i) res[[i]][[1]]))
  xi.var.list = NULL
  if(var.est){
    ## estimated xi covariance matrice; KxK for each idseval
    xi.var.list = lapply(1:length(res), function(i) res[[i]][[2]])
  }
  ## estimated time effect matrix
  smat = t(ev.vec %*% xi.est)
  svarmat = NULL
  if(var.est){
    ## estimated time effect variance matrix
    svarmat = t(sapply(1:length(xi.var.list), 
                       function(i) apply(ev.vec, 1, 
                                         function(x) sum(x*t(x*xi.var.list[[i]])))))
  }
  return(list(
    smat = smat,
    svarmat = svarmat,
    xi.est = xi.est,
    xi.var.list = xi.var.list))
}

PACE1d = function(ids, doy, resid, ev.vec, nugg, ev.val, doyeval, idseval, var.est){
  if(missing(idseval))
    idseval = sort(unique(ids))
  if(!all(idseval %in% ids))
    stop("idseval is not in ids.")
  i = NULL
  res = foreach(i = 1:length(idseval)) %dopar%{
    id.idx = which(ids == idseval[i])
    doy.idx = which(doyeval %in% doy[id.idx])
    phi = ev.vec[doy.idx,, drop=FALSE]
    H = ev.val * t(phi)
    sigma = phi %*% H + diag(nugg, length(id.idx))
    list(
      ## estimation of xi
      H %*% solve(sigma, resid[id.idx]),
      ## variance estimation of xi
      H %*% solve(sigma, t(H))
    )
  }
  ##estimated xi matrix: K x length(idseval)
  xi.est = do.call('cbind', lapply(1:length(res), function(i) res[[i]][[1]]))
  xi.var.list = NULL
  if(var.est){
    ## estimated xi covariance matrice; KxK for each idseval
    xi.var.list = lapply(1:length(res), function(i) res[[i]][[2]])
  }
  ## estimated time effect matrix
  tmat = t(ev.vec %*% xi.est)
  tvarmat = NULL
  if(var.est){
    ## estimated time effect variance matrix
    tvarmat = t(sapply(1:length(xi.var.list), 
                       function(i) apply(ev.vec, 1, 
                                         function(x) sum(x*t(x*xi.var.list[[i]])))))
  }
  
  return(list(idseval = idseval,
              tmat = tmat,
              tvarmat = tvarmat,
              xi.est = xi.est,
              xi.var.list = xi.var.list))
}
