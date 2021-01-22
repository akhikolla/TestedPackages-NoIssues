library(minerva)
library(parallel)
setClass("MIC",

  contains = "GeneralTest"

)

setValidity("MIC", function(object){
  if(object@p.opt == "dist")
    stop('No "dist" option for MIC, please use "MC" ro "table".')

})


setMethod("test", signature(object = "MIC"), function(object){

  p = object@pdata
  data = p[[ls(p)]]
  n = nrow(data)




  # Test statistic

  MIC = minerva::mine(data[,1],data[,2])$MIC

  #table
  if(object@p.opt == "table"){
    if(n > nrow(MIC_cutoff)*100 || n < 100)
      stop('Input size must be 100 to 5000 for p.opt = "table". Use exact method by setting p.opt = "MC".')
    pv = checkTable(MIC, n, "MIC")

  }

  #MC
  else{
    #estimate run time
    if(n>50){
      est = object@num.MC*0.001269 + 3
      warning(paste("Estimated computing time: ", round(est, 1), " seconds"), immediate. = TRUE)
    }

    sn = object@num.MC

    MC = function(i){
      if(object@set.seed){set.seed(i)}
      return(mine(rnorm(n), rnorm(n))$MIC)
    }

    cl <- parallel::makeCluster(10)

    parallel::clusterExport(cl, "MC", envir = environment())
    parallel::clusterExport(cl, "n", envir = environment())
    parallel::clusterEvalQ(cl, library(minerva))
    results <- parallel::parLapply(cl, 1:sn, MC)
    MICs <- as.vector(do.call('rbind',results))
    parallel::stopCluster(cl)

    NGE = sum(MICs>MIC)
    pv = NGE/(sn+1)

  }

    #BS
    if(object@BS.CI == 0)
       return(new("testforDEP_result", TS = MIC, p_value = pv))


    times = 1000

    BS = function(i){
      if(object@set.seed){set.seed(i)}
      BSdata = data[sample(1:n, n, replace = TRUE),]
      return(mine(BSdata[,1], BSdata[,2])$MIC)
    }

     cl <- parallel::makeCluster(10)

      parallel::clusterExport(cl, "data", envir = environment())
      parallel::clusterExport(cl, "n", envir = environment())
      parallel::clusterExport(cl, "BS", envir = environment())
      parallel::clusterEvalQ(cl, library(minerva))
      results <- parallel::parLapply(cl, 1:times, BS)
      MICs <- do.call('rbind',results)
      parallel::stopCluster(cl)

      CI = getCI(MIC, MICs, object@BS.CI)
      return(new("testforDEP_result", TS = MIC, p_value = pv, CI = CI))

})
