
setClass("SPEARMAN",

  contains = "GeneralTest"

)

setValidity("SPEARMAN", function(object){
  if(object@p.opt == "table")
    stop('No "table" option for SPEARMAN, please use "MC" or "dist".')

})


setMethod("test", signature(object = "SPEARMAN"), function(object){
  p = object@pdata
  n = nrow(p[[ls(p)]])
  if(length(unique(p[[ls(p)]][,1])) != n || length(unique(p[[ls(p)]][,2])) != n)
    stop("Ties detected. Spearman cannot handle tie.")

  out = cor.test(p[[ls(p)]][,1], p[[ls(p)]][,2], method = "spearman")


  #Test Statistic
  rho = out$estimate
  t = rho*sqrt((n-2)/(1-rho^2))
  #p-value
  if(object@p.opt == "dist"){
    pv = out$p.value
  }

  else{
    sn = object@num.MC
    ts = c()
    for(i in 1:sn){
      if(object@set.seed){set.seed(i)}
      sim = data.frame(cbind(rnorm(n,0,1), rnorm(n,0,1)))
      rhos = cor.test(sim[,1], sim[,2], method = "spearman")$estimate
      ts[i] = rhos*sqrt((n-2)/(1-rhos^2))
    }
    NGE = length(which(ts>t))
    NLE = sn-NGE
    pv = 2*min(NGE/(sn+1), NLE/(sn+1))
  }

  if(object@BS.CI != 0){
    warning("Spearman can't handle tie. Bootstrap CI doesn't apply.")
  }

  return(new("testforDEP_result", TS = t, p_value = pv))

})
