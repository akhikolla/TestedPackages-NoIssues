library(Hmisc)
library(parallel)
setClass("HOEFFD",

  contains = "GeneralTest"

)

setValidity("HOEFFD", function(object){
  if(object@p.opt == "dist")
    stop('No "dist" option for HOEFFD, please use "MC" ro "table".')


})


setMethod("test", signature(object = "HOEFFD"), function(object){
  p = object@pdata
  data = p[[ls(p)]]
  n = nrow(data)

  out = Hmisc::hoeffd(data[,1], data[,2])
  D = out$D[1,2]

  #p value by table
  if(object@p.opt == "table"){
    pv = out$P[1,2]
  }

  #p value by MC
  else{
      sn = object@num.MC

      MC = function(i){
        if(object@set.seed){set.seed(i)}
        return(Hmisc::hoeffd(rnorm(n), rnorm(n))$D[1,2])
      }

      cl <- parallel::makeCluster(10)

      tryCatch({
        parallel::clusterExport(cl, "MC", envir = environment())
        parallel::clusterExport(cl, "n", envir = environment())
        parallel::clusterEvalQ(cl, library(Hmisc))
        results <- parallel::parLapply(cl, 1:sn, MC)
        Ds <- as.vector(do.call('rbind',results))

        NGE = length(which(Ds>D))
        pv = NGE/(sn+1)

        parallel::stopCluster(cl)
      }, error = function(e){
                  parallel::stopCluster(cl)
                  return(e)
              }
      )
    }
  #BS
      if(object@BS.CI != 0){
        times = 1000

        BS = function(i){
          if(object@set.seed){set.seed(i)}
          BSdata = data[sample(1:n, n, replace = TRUE),]
          return(Hmisc::hoeffd(BSdata[,1],BSdata[,2])$D[1,2])
        }

         cl <- parallel::makeCluster(10)

        tryCatch({
          parallel::clusterExport(cl, "BS", envir = environment())
          parallel::clusterExport(cl, "data", envir = environment())
          parallel::clusterExport(cl, "n", envir = environment())
          parallel::clusterEvalQ(cl, library(Hmisc))
          results <- parallel::parLapply(cl, 1:times, BS)
          Ds <- do.call('rbind',results)
          parallel::stopCluster(cl)

          CI = getCI(D, Ds, object@BS.CI)
          return(new("testforDEP_result", TS = D, p_value = pv, CI = CI))

        } , error = function(e){
                    parallel::stopCluster(cl)
                    return(e)
                }
        )
       }

      else{
        return(new("testforDEP_result", TS = D, p_value = pv))
      }





})
