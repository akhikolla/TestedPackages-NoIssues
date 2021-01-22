## fgasp Class
setClass("fgasp", 		
         representation( 
           num_obs = "integer",          ## observations number
           have_noise="logical",            ## whether the data contain noise
           kernel_type="character",           ## information of the kernel
           ## data
           input="vector",              ## location of the sort input
           delta_x="vector",            ## distance between each sorted input
           output = "vector"           ## the observations, size nx1
         ),
)

## show
if(!isGeneric("show")) {
  setGeneric(name = "show",
             def = function(object) standardGeneric("show")
  )
}

setMethod("show", "fgasp",
          function(object){
            show.fgasp(object)
          }
)



## pred.fgasp Class
setClass("predictobj.fgasp", representation(
  num_testing="vector",     ##the number of testing input
  testing_input="vector",   ##sorted testing input
  param="vector",            ##param
  mean = "vector",          ##predictive mean 
  var="vector",              ##predictive variance
  var_data="logical"        ##whether to calculate the predictive variance of the data.
),
)

if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object,...) standardGeneric("predict")
  )
}

setMethod("predict", "fgasp",
          definition=function(param=NA,object, testing_input, var_data=TRUE,sigma_2=NULL){
            predict.fgasp(param=param,object=object,testing_input=testing_input,var_data=var_data,
                          sigma_2=sigma_2)
          }
)

