
##########################################################################
## rcalibration function for multiple sources data
##########################################################################
##for the design and observations are both list

rcalibration_MS <- function(design, observations, p_theta=NULL, index_theta=NULL,
                            X=as.list(rep(0,length(design))), have_trend=rep(FALSE,length(design)),
                            simul_type=rep(0, length(design)), input_simul=NULL, output_simul=NULL,
                            simul_nug=rep(FALSE,length(design)),math_model=NULL,
                            theta_range=NULL, sd_proposal_theta=rep(0.05,p_theta),
                            sd_proposal_cov_par=NULL,
                            S=10000,S_0=1000, discrepancy_type=rep('S-GaSP',length(design)),
                            kernel_type=rep('matern_5_2',length(design)),
                            tilde_lambda=rep(1/2,length(design)),  a=NULL,b=NULL,alpha=NULL,
                            output_weights=NULL){
  
  
  if( !is.list(design) | !is.list(observations) ){
    stop("The design and observations should both be a list for data from multiple sources \n")
    
  }
  
  if(is.null(p_theta)){
    stop("Please specify `p_theta', the number of parameters to be calibrated \n")
  }
  
  num_sources=length(design)
  
  
  
  
  model <- new("rcalibration_MS")
  
  model@num_sources=num_sources
  
  ##default choice is that the all calibration parameters are shared across all sources
  if(is.null(index_theta)[1]){
    model@index_theta=as.list(1:num_sources)
    for(i in 1:num_sources){
      model@index_theta[[i]]=1:p_theta
    }
  }
  
  
  if(is.null(theta_range) ){
    stop("Please specify theta_range \n")
  }
  
  model@p_x=rep(0,num_sources)
  
  model@input=design
  model@output=observations
  
  model@output_weights=as.list(rep(0,num_sources))
  model@X=as.list(1:num_sources)
  
  for(i in 1: num_sources){
    model@p_x[i]=dim(design[[i]])[2]
    model@output_weights[[i]]=rep(1, length(model@output[[i]]))
    model@X[[i]]=as.matrix(X[[i]])
    
    if(simul_type[i]==1){
      if(is.null(math_model)){
        stop("Please specify the math_model  \n")
      }
    }
    if ( (discrepancy_type[i]!='no-discrepancy')&&(discrepancy_type[i]!='GaSP')&& (discrepancy_type[i]!='S-GaSP') ){
      stop("simul_type should be chosen from no-discrepancy, GaSP or S-GaSP \n")
    }
    if(simul_type[i]==0){
      if(is.null(input_simul[i]) | is.null(output_simul[i])  ){
        stop("input_simul and output_simul are needed to be specified for constructing the emulator \n")
        
      }
    }
  }
  
  if(!is.null(output_weights)){
    model@output_weights=output_weights
  }
  
  if(is.null(b)){
    b=rep(1, num_sources)
  }
  
  
  model@num_obs=rep(0,num_sources)
  model@theta_range=theta_range
  model@p_theta=as.integer(p_theta)
  model@have_trend=have_trend
  model@q=rep(0, num_sources)
  model@kernel_type=kernel_type
  model@alpha=as.list(rep(0,num_sources))
  model@discrepancy_type=discrepancy_type
  
  model@have_trend=have_trend
  
  model@R0 =as.list(1:num_sources)
  
  model@prior_par=as.list(rep(0,num_sources))
  
  
  #model@prior_par=c(CL,a,b)
  
  for(i in 1:num_sources){
    model@num_obs[i]=length(model@output[[i]])
    
    if(model@discrepancy_type[i]!='no-discrepancy'){
      
      CL = rep(0,model@p_x[i])    ###CL is also used in the prior so I make it a model parameter
      for(i_cl in 1:model@p_x[i]){
        CL[i_cl] = (max(model@input[[i]][,i_cl])-min(model@input[[i]][,i_cl]))/model@num_obs[i]^{1/model@p_x[i]}
      }
      model@prior_par[[i]]=CL
      
      if(is.null(a)){
        model@prior_par[[i]]=c( model@prior_par[[i]],1/2-model@p_x[i])
      }else{
        model@prior_par[[i]]=c( model@prior_par[[i]],a[i])
      }
      model@prior_par[[i]]=c( model@prior_par[[i]],b[i])
    }
    
    model@num_obs[i]=dim(design[[i]])[1]
    if(model@have_trend[i]==TRUE){
      model@q[i]=dim(model@X[[i]])[2]
    }
    
    model@R0[[i]]=as.list(1:model@p_x[i])
    
    for(i_p_x in 1:model@p_x[i]){
      model@R0[[i]][[i_p_x]] = as.matrix(abs(outer(model@input[[i]][,i_p_x], model@input[[i]][,i_p_x], "-")))
    }
    
    if(kernel_type[i]=='pow_exp'){
      model@alpha[[i]]=alpha[[i]]
    }
  }
  
  
  model@sd_proposal_theta=sd_proposal_theta
  
  model@sd_proposal_cov_par=as.list(rep(0,num_sources))
  
  
  if(is.null(sd_proposal_cov_par)){
    for(i in 1:num_sources){
      model@sd_proposal_cov_par[[i]]=c(rep(0.25,dim(design[[i]])[2]),0.25)
    }
  }else{
    model@sd_proposal_cov_par=sd_proposal_cov_par
  }
  
  
  model@tilde_lambda=tilde_lambda
  model@S=as.integer(S)
  model@S_0=as.integer(S_0)
  
  
  #prior parameter for jointly robust prior
  
  
  
  model@simul_type=as.integer(simul_type)
  
  emulator=as.list(rep(0,num_sources))
  par_cur_theta=rowMeans(theta_range)  ###current par
  
  par_cur_individual=as.list(rep(0,num_sources))
  math_model_MS=list()
  
  for(i in 1:num_sources){
    if(model@simul_type[i]==1){
      math_model_MS[[i]]=math_model[[i]]
    }
    #model@simul_type[i]=as.integer(simul_type[i])
    
    ##
    if(model@simul_type[i]==0){
      cat('Starting emulation \n')
      emulator[[i]]=rgasp(design=input_simul[[i]], response=output_simul[[i]],nugget.est=simul_nug[i])
      
      cat('end emulation \n')
      cat('The inverse range  parameters from the emulation are ', emulator[[i]]@beta_hat, '\n')
      cat('The nugget  parameter from the emulation is ', emulator[[i]]@nugget, '\n')
      model@emulator=emulator
    }
    
    if(discrepancy_type[i]=='no-discrepancy'){
      par_cur_individual[[i]]= 1 #var par
    }else{
      par_cur_individual[[i]]=log(10/colMeans(abs(model@input[[i]])) )
      par_cur_individual[[i]]= c(par_cur_individual[[i]],0,1) #nugget and var par
    }
    
    if(model@have_trend[[i]]){
      par_cur_individual[[i]]=c( par_cur_individual[[i]],rep(0,dim(model@X[[i]])[2]))
    }
    
  }
  
  
  #cat('success \n')
  #cat(p_theta,' \n')
  
  ###initialize MCMC
  
  
  sample_list=post_sample_MS(model,par_cur_theta, par_cur_individual, emulator,math_model_MS)
  
  
  
  
  model@post_theta=as.matrix(sample_list$record_theta[((S_0+1):S),]) 
  #model@post_individual_par=as.list(rep(0,model@num_sources))
  model@post_individual_par=list()
  
  index_emulator=0
  
  for(i in 1:model@num_sources){
    model@post_individual_par[[i]]=sample_list$individual_par[[i]][((S_0+1):S),] 
    if(model@simul_type[i]==0){
      index_emulator=1
    }
  }
  model@post_value=sample_list$record_post[((S_0+1):S),] 
  
  ##model@accept_S=c(sample_list$accept_S_theta,sample_list$accept_S_beta)
  
  model@accept_S_theta=sample_list$accept_S_theta
  model@accept_S_beta=sample_list$accept_S_beta
  
  model@count_boundary=sum(sample_list$count_dec_record)
  
  #sample_list
  return(model)
  
  # if(index_emulator==1){## if there is an emulator, we also return it
  #   return_list=list()
  #   return_list$model=model
  #   return(return_list)
  # }else{
  #   return(model)
  # }
}
