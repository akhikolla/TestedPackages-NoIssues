
setClass("GeneralTest",

  slots = list(
    pdata = "environment",
    p.opt = "character",
    num.MC = "numeric",
    BS.CI = "numeric",
    set.seed = "logical"
    ),

  prototype = list(
    p.opt = "dist",
    num.MC = 10000,
    BS.CI = 0,
    set.seed = FALSE
  )

)


setValidity("GeneralTest", function(object){

  p = object@pdata
  if(!is.data.frame(p[[ls(p)]])){
    stop("Data must be a data frame.")
  }

  if(ncol(p[[ls(p)]]) != 2 || nrow(p[[ls(p)]][][1]) != nrow(p[[ls(p)]][2])){
    delete("p", environment())
    stop("Data must be two-columned, equal-length data frame")
  }

  names(p[[ls(p)]]) = c("x", "y")

  if(!is.numeric(p[[ls(p)]]$x) || !is.numeric(p[[ls(p)]]$y)){
    delete("p", environment())
    stop('Data must be type: "numeric".')
  }

  if(object@p.opt != "dist" && object@p.opt != "MC" && object@p.opt != "table"){
    delete("p", environment())
    stop('Please specify p-value option: "dist", "MC", or "table".')
  }

  if(object@num.MC < 100 ||object@num.MC > 10000){
    delete("p", environment())
    stop('Please specify Monte Carlo simulation numbers (100 to 10000).')
  }

  if(object@BS.CI < 0 || object@BS.CI > 1){
    delete("p", environment())
    stop("BS.CI must be numeric between 0 and 1.")
  }


})

setGeneric("test", function(object){standardGeneric("test")})
