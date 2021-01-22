


predict.rcalibration<-function(object, testing_input, X_testing=matrix(0,dim(testing_input)[1],1),
                               n_thinning=10, testing_output_weights=rep(1,dim(testing_input)[1]), 
                               interval_est=NULL,interval_data=F,math_model=NULL,...){
  predictobj<- new("predictobj.rcalibration")
  record_cm_pred=0
  record_cm_pred_no_trend=0
  
  if(object@discrepancy_type=='no-discrepancy'){
    if(!is.null(interval_est)){
      record_interval=matrix(0,dim(testing_input)[1],length(interval_est));
    }
    SS=floor(dim(object@post_sample)[1]/n_thinning)
    c_prop=1/4
    
    for(i_S in (1: SS)*n_thinning ){
      #print(i_S)
      
        if(i_S==floor(SS*n_thinning*c_prop)){
          cat(c_prop*100, 'percent is completed \n')
          c_prop=c_prop+1/4
        }
        
     # print(i_S)
      
      theta=object@post_sample[i_S,1:object@p_theta]
      sigma_2_delta=object@post_sample[i_S,object@p_theta+1]

      mean_cm_test=mathematical_model_eval(testing_input,theta,object@simul_type,object@emulator,math_model);
      mean_cm_test_no_trend=mean_cm_test
      if(object@have_trend){
        theta_m=object@post_sample[i_S,(object@p_theta+2):(object@p_theta+1+object@q)]
        mean_cm_test=mean_cm_test+X_testing%*%theta_m
      }
      
      
      record_cm_pred=record_cm_pred+mean_cm_test
      
      record_cm_pred_no_trend=record_cm_pred_no_trend+mean_cm_test_no_trend
      if(!is.null(interval_est)){
        qnorm_all=qnorm(interval_est);
        for(i_int in 1:length(interval_est) ){
          record_interval[,i_int]=record_interval[,i_int]+mean_cm_test+qnorm_all[i_int]*sqrt(sigma_2_delta/testing_output_weights)
        }
      }
      
    }
    record_cm_pred=record_cm_pred/SS
    record_cm_pred_no_trend=record_cm_pred_no_trend/SS
    #output.list <- list()
  
    #output.list$math_model_mean=record_cm_pred
    
    predictobj@math_model_mean=record_cm_pred
    predictobj@math_model_mean_no_trend=record_cm_pred_no_trend
    
    if(!is.null(interval_est)){
      record_interval=record_interval/SS
      #output.list$interval=record_interval
      predictobj@interval=record_interval
      #ans.list[[2]]=record_interval
    }
    
    return(predictobj)
    
  }else{
  if(!is.null(interval_est)){
    c_star_record=rep(0,dim(testing_input)[1])
    record_interval=matrix(0,dim(testing_input)[1],length(interval_est));
  }
  
  N_testing=dim(testing_input)[1]
  r0=as.list(1:object@p_x)
  for(i in 1:object@p_x){
    r0[[i]]=abs(outer(object@input[,i],testing_input[,i],'-'))
  }
  
  record_pred_mean=0
  
  SS=floor(dim(object@post_sample)[1]/n_thinning)
  
  #c_star=rep(0, n_testing)
  c_prop=1/4
  
  for(i_S in (1: SS)*n_thinning ){
    if(i_S==floor(SS*n_thinning*c_prop)){
      cat(c_prop*100, 'percent is completed \n')
      c_prop=c_prop+1/4
    }
    
    theta=object@post_sample[i_S,1:object@p_theta]
    beta_delta=exp(object@post_sample[i_S,(object@p_theta+1):(object@p_theta+object@p_x)])
    eta_delta=exp(object@post_sample[i_S,object@p_theta+object@p_x+1])
    sigma_2_delta=object@post_sample[i_S,object@p_theta+object@p_x+2]
    if(object@have_trend){
      theta_m=object@post_sample[i_S,(object@p_theta+object@p_x+3):(object@p_theta+object@p_x+2+object@q)]
    }
    
    mean_cm=mathematical_model_eval(object@input,theta,object@simul_type,object@emulator,math_model);

    output_minus_cm=object@output- mean_cm
    
    if(object@have_trend){
       output_minus_cm=output_minus_cm-object@X%*%theta_m
    }
    ##remember 
    if(object@discrepancy_type=="GaSP"){
      L=Get_R_new( (beta_delta),   (eta_delta),  
                   object@R0, object@kernel_type,object@alpha,1/object@output_weights)
    }else if(object@discrepancy_type=="S-GaSP"){
      L=Get_R_z_new( (beta_delta),   (eta_delta),  object@tilde_lambda,
                     object@R0, object@kernel_type,object@alpha,1/object@output_weights)
    }

    R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm ))
    
    if(object@discrepancy_type=="S-GaSP"){
      R_inv_y=Update_R_inv_y(R_inv_y,object@R0,beta_delta,object@kernel_type,object@alpha,object@tilde_lambda,object@num_obs)
    }
    
    r=separable_kernel(r0,beta_delta,kernel_type=object@kernel_type,object@alpha)
    
    rt_R_inv_y=t(r)%*%R_inv_y
    
    mean_cm_test=mathematical_model_eval(testing_input,theta,object@simul_type,object@emulator,math_model);
    
    mean_cm_test_no_trend=mean_cm_test
      #f_M_testing=mean_cm_test
    
    if(object@have_trend){
      mean_cm_test=mean_cm_test+X_testing%*%theta_m
    }
    
    
    pred_mean=mean_cm_test+rt_R_inv_y
    
    record_cm_pred=record_cm_pred+mean_cm_test
    record_cm_pred_no_trend=record_cm_pred_no_trend+mean_cm_test_no_trend
    
    record_pred_mean=record_pred_mean+pred_mean
    
    if(!is.null(interval_est)){
      if(object@discrepancy_type=="GaSP"){
        R_tilde_inv_r=(backsolve(t(L),forwardsolve(L, r )))
        
        for(i_testing in 1:N_testing){
          c_star_record[i_testing]=1-r[,i_testing]%*%R_tilde_inv_r[,i_testing]
        }
        var_gasp_f=sigma_2_delta*c_star_record
        if(interval_data==T){
          var_gasp_f=var_gasp_f+sigma_2_delta*eta_delta/testing_output_weights
        }
          
        qnorm_all=qnorm(interval_est);
        for(i_int in 1:length(interval_est) ){
          record_interval[,i_int]=record_interval[,i_int]+pred_mean+qnorm_all[i_int]*sqrt(var_gasp_f)
          
        }
        
      }else if(object@discrepancy_type=="S-GaSP"){
        R=separable_kernel(object@R0, beta_delta, object@kernel_type,object@alpha)
        R_tilde=R+1/object@tilde_lambda*diag(object@num_obs);


        #system.time(Chol_Eigen(R_tilde))
        
        L_R_tilde=Chol_Eigen(R_tilde)
        
        I_minus_R_R_tilde=diag(object@num_obs)-t(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde,R )))
        
        
        R_tilde_inv_r=(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde, r )))

        r_z=I_minus_R_R_tilde%*%r
        
        R_z_tilde_inv_r_z=(backsolve(t(L),forwardsolve(L, r_z )))
        
        for(i_testing in 1:N_testing){
          c_star_record[i_testing]=1-r[,i_testing]%*%R_tilde_inv_r[,i_testing]-r_z[,i_testing]%*%R_z_tilde_inv_r_z[,i_testing]
        }
        var_gasp_f=sigma_2_delta*c_star_record
        
        if(interval_data==T){
          var_gasp_f=var_gasp_f+sigma_2_delta*eta_delta/testing_output_weights
        }
        
        qnorm_all=qnorm(interval_est);
        for(i_int in 1:length(interval_est) ){
          record_interval[,i_int]=record_interval[,i_int]+pred_mean+qnorm_all[i_int]*sqrt(var_gasp_f)
          
        }
        
      }
      
    }
    
    
    
  }
  
  record_cm_pred=record_cm_pred/SS
  record_pred_mean=record_pred_mean/SS
  record_cm_pred_no_trend=record_cm_pred_no_trend/SS
  
  #output.list <- list()
  
  #output.list$math_model_mean=record_cm_pred
  #output.list$mean=record_pred_mean
  
  predictobj@math_model_mean=record_cm_pred
  predictobj@math_model_mean_no_trend=record_cm_pred_no_trend
  predictobj@mean=record_pred_mean
  #ans.list=as.list(1:3)
  
  #ans.list[[1]]=record_cm_pred
  #ans.list[[2]]=record_pred_mean
  #ans.list[[3]]=NULL
  if(!is.null(interval_est)){
    record_interval=record_interval/SS
    predictobj@interval=record_interval
  }

  
  return(predictobj)
  }
}

  


predict_discrepancy_separable_2dim<-function(object, testing_input_separable,
X_testing=matrix(0,length(testing_input_separable[[1]])*length(testing_input_separable[[2]]),1 ),
n_thinning=10,interval_est=NULL,math_model=NULL,...){

  predictobj<- new("predictobj.rcalibration")
  
  r0_separable=as.list(1:object@p_x)

  for(i_x in 1:object@p_x){
     r0_separable[[i_x]]=abs(outer(object@input[,i_x],testing_input_separable[[i_x]],'-'))
  }
  
  testing_input=expand.grid(testing_input_separable[[1]],testing_input_separable[[2]])
  testing_input=as.matrix(testing_input)
  
  record_cm_pred=0
  record_cm_pred_no_trend=0
  record_pred_mean=0
  
  SS=floor(dim(object@post_sample)[1]/n_thinning)
  
  r_separable=as.list(1:object@p_x)
  

    c_prop=1/4
    
  for(i_S in (1: SS)*n_thinning ){
    
    if(i_S==floor(SS*n_thinning*c_prop)){
      cat(c_prop*100, 'percent is completed \n')
      c_prop=1/4
      
    }
    
    
    print(i_S)

    theta=object@post_sample[i_S,1:object@p_theta]
    beta_delta=exp(object@post_sample[i_S,(object@p_theta+1):(object@p_theta+object@p_x)])
    eta_delta=exp(object@post_sample[i_S,object@p_theta+object@p_x+1])
    sigma_2_delta=object@post_sample[i_S,object@p_theta+object@p_x+2]

    mean_cm=mathematical_model_eval(object@input,theta,object@simul_type,object@emulator,math_model);
    
    # if(object@simul_type==2 | object@simul_type==3){
    #    mean_cm=Mogihammer(theta,input,object@simul_type)
    # }
    output_minus_cm=object@output- mean_cm
    if(object@have_trend){
      theta_m=object@post_sample[i_S,(object@p_theta+object@p_x+3):(object@p_theta+object@p_x+2+object@q)]
      output_minus_cm=output_minus_cm-object@X%*%theta_m
    }
    ##remember 
    if(object@discrepancy_type=="GaSP"){
        L=Get_R_new( (beta_delta),   (eta_delta),  
                     object@R0, object@kernel_type,object@alpha,1/object@output_weights)
    }else if(object@discrepancy_type=="S-GaSP"){
        L=Get_R_z_new( (beta_delta),   (eta_delta),  object@tilde_lambda,
                     object@R0, object@kernel_type,object@alpha,1/object@output_weights)
    }
    

      
    
    R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm ))
    
    
    
      if(object@discrepancy_type=="S-GaSP"){
        R_inv_y=Update_R_inv_y(R_inv_y,object@R0,beta_delta,object@kernel_type,object@alpha,object@tilde_lambda,object@num_obs)
      }


    if(object@kernel_type=="matern_5_2"){
      for(i_x in 1:object@p_x){
        r_separable[[i_x]]=matern_5_2_funct(r0_separable[[i_x]],beta_delta[i_x])
      }
    }else if(object@kernel_type=="matern_3_2"){
      for(i_x in 1:object@p_x){
        r_separable[[i_x]]=matern_3_2_funct(r0_separable[[i_x]],beta_delta[i_x])
      }
    }else if(object@kernel_type=="pow_exp"){
      for(i_x in 1:object@p_x){
         r_separable[[i_x]]=pow_exp_funct(r0_separable[[i_x]],beta_delta[i_x],object@alpha[i_x])
      }
    }
    
    ###so this is only for for image or p_x=2
    
    r2_dot_R_inv_y=r_separable[[2]]*as.vector(R_inv_y)
    
    rt_R_inv_y=as.vector(t(r_separable[[1]])%*%r2_dot_R_inv_y)
    
    mean_cm_test=mathematical_model_eval(testing_input,theta,object@simul_type,object@emulator,math_model);
    
    mean_cm_test_no_trend=mean_cm_test
    
    #if(object@simul_type==2 | object@simul_type==3){
    #   mean_cm_test=Mogihammer(theta,testing_input,object@simul_type)
    #}
    #f_M_testing=mean_cm_test
    
    if(object@have_trend){
      mean_cm_test=mean_cm_test+X_testing%*%theta_m
    }
    
    
    pred_mean=mean_cm_test+rt_R_inv_y
    
    record_cm_pred=record_cm_pred+mean_cm_test
    record_cm_pred_no_trend=record_cm_pred_no_trend+mean_cm_test_no_trend
    
    record_pred_mean=record_pred_mean+pred_mean
    
    
    
  }
  


  
  record_cm_pred=record_cm_pred/SS
  record_cm_pred_no_trend=record_cm_pred_no_trend/SS
  record_pred_mean=record_pred_mean/SS
  
  predictobj@math_model_mean=record_cm_pred
  predictobj@math_model_mean_no_trend=record_cm_pred_no_trend
  predictobj@mean=record_pred_mean
  #ans.list=as.list(1:3)
  
  #ans.list[[1]]=record_cm_pred
  #ans.list[[2]]=record_pred_mean
  #ans.list[[3]]=NULL
  if(!is.null(interval_est)){
    record_interval=record_interval/SS
    predictobj@interval=record_interval
  }
  
  return(predictobj)
}


