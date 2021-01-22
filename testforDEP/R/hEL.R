setClass("EL",

  contains = "GeneralTest"

)

setValidity("EL", function(object){
  if(object@p.opt == "dist")
    stop('No "dist" option for EL, please use "MC" or "table".')

})

setMethod("test", signature(object = "EL"), function(object){
  p = object@pdata
  data = p[[ls(p)]]
  n = nrow(data)




  if(length(unique(data[,1])) != n || length(unique(data[,2])) != n)
    warning("Tie detected in input. Rank (order statistics) will be computed by First occurrence.", immediate. = TRUE)


  #test statistics
  tn = Tn(data[,1], data[,2])

  #table
  if(object@p.opt == "table"){
    if(n > nrow(EL_cutoff)*100 || n < 100)
      stop('Input size must be 100 to 5000 for p.opt = "table". Use exact method by setting p.opt = "MC".')
    pv = checkTable(tn, n, "EL")

  }

  #MC
  else {

    #estimate run time
    if(n>100 && n<=500){
      est = EL_time[round(n / 100), round(object@num.MC / 1000)]
      warning(paste("Estimated computing time: ", round(est, 1), " seconds"), immediate. = TRUE)
    }
    else if(n>500)
       warning("Estimated computing time: over 3 min", immediate. = TRUE)

    MC = function(i){
      if(object@set.seed){set.seed(i)}
      return (MC_EL_count(tn, n, sn/8))
    }

    sn = object@num.MC
    cl <- parallel::makeCluster(8)

    parallel::clusterExport(cl, "tn", envir = environment())
    parallel::clusterExport(cl, "n", envir = environment())
    parallel::clusterExport(cl, "sn", envir = environment())
    parallel::clusterExport(cl, "MC", envir = environment())
    results <- parallel::parLapply(cl, 1:8, MC)
    counts <- as.vector(do.call('rbind',results))
    parallel::stopCluster(cl)

    pv = sum(counts)/(sn+1)

  }

  if(object@BS.CI == 0){
      return(new("testforDEP_result", TS = tn, p_value = pv))
  }

  #BS.CI
  times = 1000
  BS = function(i){
    if(object@set.seed){set.seed(i)}
    BSdata = data[sample(1:n, n, replace = T),]
    return (Tn(BSdata[,1], BSdata[,2]))
  }

  cl <- parallel::makeCluster(8)

  parallel::clusterExport(cl, "n", envir = environment())
  parallel::clusterExport(cl, "BS", envir = environment())
  parallel::clusterExport(cl, "data", envir = environment())
  results <- parallel::parLapply(cl, 1:times, BS)
  Tns <- as.vector(do.call('rbind',results))
  parallel::stopCluster(cl)

  CI = getCI(tn, Tns, object@BS.CI)
  return(new("testforDEP_result", TS = tn, p_value = pv, CI = CI))



})


