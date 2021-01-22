setClass("testforDEP_result",

  slots = list(
    TS = "numeric",
    p_value = "numeric",
    CI = "list"
    ),

  prototype = list(
    CI = list()
  )

)

setValidity("testforDEP_result", function(object){
  if(object@p_value >1 || object@p_value <0){
    stop("p_value is outside 0~1")
  }

})
