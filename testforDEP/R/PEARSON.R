
setClass("PEARSON",

  contains = "GeneralTest"

)

setValidity("PEARSON", function(object){
  if(object@p.opt == "table")
    stop('No "table" option for PEARSON, please use "MC" or "dist".')
})


setMethod("test", signature(object = "PEARSON"), function(object){
  p = object@pdata
  r = cor(p[[ls(p)]])[1,2]
  n = nrow(p[[ls(p)]])
  t = r*sqrt((n-2)/(1-r^2))

  #dist
  if(object@p.opt == "dist"){
    pv = 2*min(pt(t, n-2), 1-pt(t, n-2))
  }

  #MC
  else{
    sn = object@num.MC
    ts = c()
    for(i in 1:sn){
      if(object@set.seed){set.seed(i)}
      sim = cbind(rnorm(n,0,1), rnorm(n,0,1))
      r = cor(sim)[1,2]
      ts[i] = r*sqrt((n-2)/(1-r^2))
    }

    NGE = length(which(ts>t))
    NLE = sn-NGE
    pv = 2*min(NGE/(sn+1), NLE/(sn+1))

  }

  #BS
  if(object@BS.CI != 0){
    times = 1000
    ts = c()
    for(i in 1:times){
      if(object@set.seed){set.seed(i)}
      index = sample(1:n, n, replace = TRUE)
      r = cor(p[[ls(p)]][index,])[1,2]
      ts[i] = r*sqrt((n-2)/(1-r^2))
    }

    CI = getCI(t, ts, object@BS.CI)
    return(new("testforDEP_result", TS = t, p_value = pv, CI = CI))
  }

  else{
    return(new("testforDEP_result", TS = t, p_value = pv))
  }
})


