
setClass("KENDALL",

  contains = "GeneralTest"

)

setValidity("KENDALL", function(object){
  if(object@p.opt == "table")
    stop('No "table" option for KENDALL, please use "MC" or "dist".')
})


setMethod("test", signature(object = "KENDALL"), function(object){
  p = object@pdata
  n = nrow(p[[ls(p)]])
  out = cor.test(p[[ls(p)]][,1], p[[ls(p)]][,2], method = "kendall")
  #Test Statistic
  z = as.numeric(out$statistic)
  #p-value
  if(object@p.opt == "dist"){
    pv = out$p.value
  }

  else{
    sn = object@num.MC
    zs = c()
    for(i in 1:sn){
      if(object@set.seed){set.seed(i)}
      sim = data.frame(cbind(rnorm(n,0,1), rnorm(n,0,1)))
      zs[i] = cor.test(sim[,1], sim[,2], method = "kendall")$statistic
    }
    NGE = length(which(zs>z))
    NLE = sn-NGE
    pv = 2*min(NGE/(sn+1), NLE/(sn+1))
  }

  if(object@BS.CI != 0){
    times = 1000
    zs = c()
    for(i in 1:times){
      if(object@set.seed){set.seed(i)}
      index = sample(1:n, n, replace = TRUE)
      zs[i] = cor.test(p[[ls(p)]][index,1], p[[ls(p)]][index,2], method = "kendall")$statistic
    }

    CI = getCI(z, zs, object@BS.CI)
    return(new("testforDEP_result", TS = z, p_value = pv, CI = CI))

  }

  else{
    return(new("testforDEP_result", TS = z, p_value = pv))
  }



})
