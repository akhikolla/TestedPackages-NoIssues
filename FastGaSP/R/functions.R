##########################################################################
## fgasp construction function
## 
## FastGaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2018-present Mengyang Gu
##							  
##    
##########################################################################


fgasp <- function(input, output, have_noise=TRUE, kernel_type='matern_5_2'){


  if( (kernel_type!='matern_5_2')&(kernel_type!='exp') ){
    stop("the current version only support the Matern covariance with smoothness 
         parameter being 2.5 or 0.5 (exponential kernel). \n")
  }
  
  ##if the input is not sorted so we will sort them
  object <- new("fgasp")
  object@num_obs=length(output)
  #object@param=param
  object@kernel_type=kernel_type
  object@have_noise=have_noise
  
  ## if the input is not sorted, then I sort it
  if(sum(which(input[2:object@num_obs]-input[1:(object@num_obs-1)]<0 ))>0){
    input_sorted_all=sort(input,index.return=TRUE)
    object@input=input_sorted_all[[1]]
    object@output=output[input_sorted_all[[2]]]
  }else{
    object@input=input
    object@output=output
  }
  object@delta_x=object@input[2:object@num_obs]-object@input[1:(object@num_obs-1)]
  
  if(length(object@delta_x)!=(length(object@output)-1)){
    stop("length(output) should be equal to length(delta_x)+1 \n")
  }
  ##the current version does not support duplicated inputs
  if(length(which(object@delta_x==0))>0){
    stop("Please delete the repeated inputs. \n")
  }
  return(object)
}

show.fgasp <- function(object) {	
  cat('Number of observations: ',object@num_obs,'\n')
  cat('Have noise: ',object@have_noise,'\n')
  cat('Type of kernel: ',object@kernel_type,'\n')
  
}


##likelihood function of GaSP using fast computation
##the Filtering algorithm is implemented using Rcpp and RcppEigen
log_lik<-function(param,object){
  param=as.vector(param)
  if(object@have_noise==T && length(param)!=2){
       stop("Please specify both the log inverse range parameter 
            and log nugget parameter if there is noise. \n")
  }
  
  log_det_S2=Get_log_det_S2(param,object@have_noise,object@delta_x,
                            object@output, object@kernel_type);
  ##log likelihood
  -log_det_S2[[1]]/2-(object@num_obs)/2*log(log_det_S2[[2]])
}

##prediction 
predict.fgasp<-function(param,object,testing_input,var_data=TRUE,sigma_2=NULL){
  param=as.vector(param)
  if(object@have_noise==T && length(param)!=2){
    stop("Please specify both the log inverse range parameter 
          and log nugget parameter. \n")
  }
  ## if sigma_2 is not given, I will estimate it
  if(is.null(sigma_2)){
  sigma_2=Get_log_det_S2(param,object@have_noise,object@delta_x,
                             object@output,object@kernel_type)[[2]]/object@num_obs
  }
  
  predictobj<- new("predictobj.fgasp")
  
  testing_input_sorted=sort(testing_input)
  predictobj@num_testing=length(testing_input_sorted)
  predictobj@testing_input=testing_input_sorted
  predictobj@param=param
  predictobj@var_data=var_data
  ## we don't support the repeated testing input in this version
  if(length(which(testing_input_sorted[2:predictobj@num_testing]-
                  testing_input_sorted[1:(predictobj@num_testing-1)]==0))>0){
    stop("Please delete the repeated testing inputs. \n")
  }
  ##combine the input and testing input
  input_all_ori=c(object@input,testing_input_sorted)
  input_all_sort=sort(input_all_ori,index.return=TRUE)

  input_all=input_all_sort[[1]] ##all sorted input with and without the observations
  delta_x_all=input_all[2:length(input_all_ori)]-input_all[1:(length(input_all_ori)-1)]
  
  ##create a vector of inputs where the one mean there is an observation
  ##and zero means no observation
  index_all=rep(0, length(input_all))
  index_all[1:object@num_obs]=1
  index_all=index_all[input_all_sort[[2]]]
  index_all=as.integer(index_all)
  
  

  testing_loc=which(index_all==0)
  delta_x_all_here=delta_x_all
  
  if(length(which(delta_x_all==0))>0){
    #input_all=input_all[-which(delta_x_all==0)]
    index_all=index_all[-(which(delta_x_all==0)+1)]
    delta_x_all_here=delta_x_all[-which(delta_x_all==0)]
  }
  
  
  
  KF_smoother_result=Kalman_smoother(param, object@have_noise,
                                     index_all, delta_x_all_here,object@output,sigma_2,object@kernel_type)
  if(length(which(delta_x_all==0))==0){
    predictobj@mean=KF_smoother_result[[1]][testing_loc]
    if(var_data==T | object@have_noise==F){
       predictobj@var=KF_smoother_result[[2]][testing_loc]
    }else{
      predictobj@var=KF_smoother_result[[2]][testing_loc]-sigma_2*exp(param[2])
    }
  }else{##enlarge it to the orginal long sequence that may contain knots. 
    res=rep(NA,length(input_all_ori))
    res[-(which(delta_x_all==0)+1)]=KF_smoother_result[[1]]
    res[which(delta_x_all==0)+1]= res[which(delta_x_all==0)]
    predictobj@mean=res[testing_loc]
    
    res[-(which(delta_x_all==0)+1)]=KF_smoother_result[[2]]
    res[which(delta_x_all==0)+1]= res[which(delta_x_all==0)]
    if(var_data==T | object@have_noise==F){
      predictobj@var=res[testing_loc]
    }else{
      predictobj@var=res[testing_loc]-sigma_2*exp(param[2])
    }
    
  }
  
  return(predictobj)
  
}

