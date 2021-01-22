##########################################################################
## rcalibration fit function
##########################################################################


rcalibration <- function(design, observations, p_theta=NULL,X=matrix(0,dim(design)[1],1), have_trend=FALSE,
                         simul_type=0, input_simul=NULL, output_simul=NULL,simul_nug=FALSE,math_model=NULL,
                         theta_range=NULL, sd_proposal=c(rep(0.05,p_theta),rep(0.25,dim(design)[2]),0.25),
                         S=10000,S_0=1000, discrepancy_type='S-GaSP', kernel_type='matern_5_2',
                         tilde_lambda=1/2,  a=1/2-dim(design)[2],b=1,alpha=rep(1.9,dim(design)[2]), output_weights=rep(1,dim(design)[1])){
  
  p_theta=as.integer(p_theta)
  
  if(simul_type==1){
    if(is.null(math_model) | is.null(theta_range)){
      stop("Please specify the math_model and theta_range \n")
    }
  }
  if(!is.null(math_model) ){ ##in this case the math model is given
    simul_type=1;
    if(is.null(theta_range) ){
      stop("Please specify the theta_range \n")
      
    }
  }
  if(is.null(p_theta)){
    stop("Please specify `p_theta', the number of parameters to be calibrated \n")
  }
  
  if (simul_type>3 | simul_type<0){
    stop("simul_type should be chosen from 0, 1, 2 or 3\n")
  }
  if ( (discrepancy_type!='no-discrepancy')&&(discrepancy_type!='GaSP')&& (discrepancy_type!='S-GaSP') ){
    stop("simul_type should be chosen from no-discrepancy, GaSP or S-GaSP \n")
  }
  
  if(simul_type==0){
    if(is.null(input_simul) | is.null(output_simul)  ){
      stop("input_simul and output_simul are needed to be specified for constructing the emulator \n")
      
    }
  }
  
  
  
  
  model <- new("rcalibration")
  model@p_x=dim(design)[2]
  #model@p_theta=p_theta
  model@num_obs=dim(design)[1]
  
  if(simul_type==0){
    model@p_theta=dim(input_simul)[2]-model@p_x
    if(model@p_theta!=p_theta){
      stop("the dimension of the input_simul is wrong \n");
    }
    ##theta_range=matrix(0,model@p_theta,2);
    if(is.null(theta_range)){
      theta_range[,1]=apply(input_simul[(model@p_x+1):(model@p_x+model@p_theta),],2,min)
      theta_range[,1]=theta_range[,1]-0.001*abs(theta_range[,1]) ##do not let the design point be the boundary
      theta_range[,2]=apply(input_simul[(model@p_x+1):(model@p_x+model@p_theta),],2,max)
      theta_range[,2]=theta_range[,2]+0.001*abs(theta_range[,2])
    }
    
  }else if (simul_type==1){
    model@p_theta=p_theta;
  }else if(simul_type==2 | simul_type==3){ ##Kilauea volcano
    model@p_theta=as.integer(5);
    if(is.null(theta_range)){
      theta_range=matrix(0,model@p_theta,2)
      theta_range[1,]=c(-2000, 3000)  
      theta_range[2,]=c(-2000, 5000)  
      theta_range[3,]=c(500, 6000)  
      theta_range[4,]=c(0, .3)   
      theta_range[5,]=c(.25, .33)  
    }
  }
  
  
  model@input=design
  model@output=observations
  model@X=X
  model@q=dim(model@X)[2]
  
  model@R0 = as.list(1:model@p_x)
  for(i in 1:model@p_x){
    model@R0[[i]] = as.matrix(abs(outer(model@input[,i], model@input[,i], "-")))
  }
  model@kernel_type=kernel_type
  model@alpha=alpha
  
  model@tilde_lambda=tilde_lambda
  model@S=as.integer(S)
  model@S_0=as.integer(S_0)
  #prior parameter for jointly robust prior
  CL = rep(0,model@p_x)    ###CL is also used in the prior so I make it a model parameter
  
  for(i_cl in 1:model@p_x){
    CL[i_cl] = (max(model@input[,i_cl])-min(model@input[,i_cl]))/model@num_obs^{1/model@p_x}
  }
  model@prior_par=c(CL,a,b)
  
  model@output_weights=output_weights
  
  
  
  model@have_trend=have_trend
  
  if(have_trend==F){
    model@X=matrix(rep(0,model@num_obs),model@num_obs,1)
    model@q=as.integer(0);
  }
  
  
  model@theta_range=theta_range
  
  model@sd_proposal=sd_proposal
  
  if(discrepancy_type=='no-discrepancy'){
    model@sd_proposal=sd_proposal[1:p_theta]
  }
  model@discrepancy_type=discrepancy_type
  
  model@simul_type=as.integer(simul_type)
  
  ##
  if(model@simul_type==0){
    cat('Starting emulation \n')
    emulator=rgasp(design=input_simul, response=output_simul,nugget.est=simul_nug)
    
    param_emulator=c(emulator@beta_hat)
    if(simul_nug==T){
      param_emulator=c(param_emulator,emulator@nugget)
    }
    model@emulator=emulator
    cat('end emulation \n')
    cat('The inverse range  parameters from the emulation are ', emulator@beta_hat, '\n')
    cat('The nugget  parameter from the emulation is ', emulator@nugget, '\n')
    
  }
  
  ###initial parameter
  if(simul_type==0){
    if(model@p_theta>1){
      #par_cur=c(colMeans(input_simul[(model@p_x+1):(model@p_x+model@p_theta),]))
      par_cur=c(colMeans(input_simul[,(model@p_x+1):(model@p_x+model@p_theta)]))
      
    }else{
      #par_cur=c(mean(input_simul[(model@p_x+1):(model@p_x+model@p_theta),]))
      par_cur=c(mean(input_simul[,(model@p_x+1):(model@p_x+model@p_theta)]))
      
    }
  }else if(simul_type==1){
    par_cur=rowMeans(theta_range)  ###current par
  }else if(simul_type==2 | simul_type==3){
    par_cur=rowMeans(theta_range)  ###current par
  }
  
  if(discrepancy_type=='no-discrepancy'){
    par_cur= c(par_cur,1) #var par
  }else{
    for(i_x in 1: model@p_x){
      par_cur= c(par_cur,log(10/sd(model@input[,i_x])))
      #   Input_list[[8]]= c(Input_list[[8]], -5)
    }
    par_cur= c(par_cur,0,1) #nugget and var par
  }
  if(model@have_trend){
    par_cur=c(par_cur,rep(0,dim(model@X)[2]))
  }
  #cat('success \n')
  #cat(p_theta,' \n')
  
  
  
  ###initialize MCMC
  sample_list=post_sample(model@input, model@output, model@R0, model@kernel_type, model@p_theta, model@output_weights,
                          par_cur, model@tilde_lambda, model@prior_par,model@theta_range,model@S, model@X, model@have_trend, model@alpha,model@sd_proposal,
                          model@discrepancy_type, model@simul_type,emulator,math_model);
  
  model@post_sample=sample_list[[1]][((S_0+1):S),]
  model@post_value=sample_list[[2]][((S_0+1):S)]
  model@accept_S=sample_list[[3]]
  model@count_boundary=sum(sample_list[[4]])
  
  return(model)
  
  # if(simul_type==0){## if there is an emulator, we also return it
  #   return_list=list()
  #   return_list$model=model
  #   return_list$emulator=emulator
  #   return(return_list)
  # }else{
  #   return(model)
  # }
  
}


show.rcalibration <- function(object) {	
  #cat("\n")
  #cat("Call:\n")
  #print(object@call)
  
  
  #cat('The dimension  of the design is: ',dim(object@input),'\n')
  #cat('The number of the output is: ',dim(object@output),'\n')
  cat('type of discrepancy function: ',object@discrepancy_type,'\n')
  cat('number of after burn-in posterior: ',object@S-object@S_0,'\n')
  cat(object@accept_S[1]/object@S, 'of proposed calibration parameters are accepted \n')
  #cat(object@accept_S[2]/object@S, 'of proposed range and nugget parameters are accepted \n')
  for( i in 1: object@p_theta){
    cat('median and 95% posterior credible interval of calibration parameter ',i, 'is', quantile(object@post_sample[,i],c(0.025,0.5,0.975) ),'\n')
  }
  
  
  # cat('The slots of this object are: ',slotNames(object),'\n')
  
}
