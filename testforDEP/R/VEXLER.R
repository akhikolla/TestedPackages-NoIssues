
library(parallel)
setClass("VEXLER",

  contains = "GeneralTest"

)

setValidity("VEXLER", function(object){
  if(object@p.opt == "dist")
    stop('No "dist" option for VEXLER, please use "MC" or "table".')

})


setMethod("test", signature(object = "VEXLER"), function(object){

  p = object@pdata
  data = p[[ls(p)]]
  n = nrow(data)

  if(length(unique(data[,1])) != n || length(unique(data[,2])) != n)
    warning("Tie detected in input. Rank (order statistics) will be computed by First occurrence.", immediate. = TRUE)





  #test statistic (log(VT))
  VT = vex(data[,1], data[,2])

  #table
  if(object@p.opt == "table"){
    if(n > nrow(vexler_cutoff)*100 || n < 100)
      stop('Input size must be 100 to 5000 for p.opt = "table". Use exact method by setting p.opt = "MC".')
    pv = checkTable(VT, n, "vexler")

  }

  #MC
  else{
    #estimate run time
    if(n>100 && n<=1000){
      est = vexler_time[round(n / 100), round(object@num.MC / 1000)]
      warning(paste("Estimated computing time: ", round(est, 1), " seconds"), immediate. = TRUE)
    }
    else if(n>1000)
       warning("Estimated computing time: over 2 min", immediate. = TRUE)


    MC = function(i){
      if(object@set.seed){set.seed(i)}
      return (MC_count(VT, n, sn/8))
    }

    sn = object@num.MC
    cl <- parallel::makeCluster(8)

    parallel::clusterExport(cl, "VT", envir = environment())
    parallel::clusterExport(cl, "n", envir = environment())
    parallel::clusterExport(cl, "sn", envir = environment())
    parallel::clusterExport(cl, "MC", envir = environment())
    results <- parallel::parLapply(cl, 1:8, MC)
    counts <- as.vector(do.call('rbind',results))
    parallel::stopCluster(cl)

    pv = sum(counts)/(sn+1)

  }


  #BS.CI

  if(object@BS.CI == 0){
        return(new("testforDEP_result", TS = VT, p_value = pv))
  }

  times = 1000
  BS = function(i){
    if(object@set.seed){set.seed(i)}
    BSdata = data[sample(1:n, n, replace = T),]
    return (vex(BSdata[,1], BSdata[,2]))
  }

  cl <- parallel::makeCluster(8)

  parallel::clusterExport(cl, "n", envir = environment())
  parallel::clusterExport(cl, "BS", envir = environment())
  parallel::clusterExport(cl, "data", envir = environment())
  results <- parallel::parLapply(cl, 1:times, BS)
  VTs <- as.vector(do.call('rbind',results))
  parallel::stopCluster(cl)

  CI = getCI(VT, VTs, object@BS.CI)
  return(new("testforDEP_result", TS = VT, p_value = pv, CI = CI))

})

