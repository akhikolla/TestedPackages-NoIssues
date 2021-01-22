
##multiple sources
predict_MS.rcalibration_MS<-function(object, testing_input, X_testing=as.list(rep(0,object@num_sources)),
                                  n_thinning=10, testing_output_weights=NULL, 
                                  interval_est=NULL,interval_data=rep(F,length(testing_input)),math_model=NULL,...){
  
  if(length(testing_input)!=object@num_sources){
    stop("The number of sources in the testing input should match the number of sources in the object. \n")
  }
  
  
  #output.list <- as.list( 1:object@num_sources)
  predictobj <- new("predictobj.rcalibration_MS")
  
  if(is.null(testing_output_weights)){
    testing_output_weights=as.list(1:2)
    for(i_source in 1: object@num_sources){
      testing_output_weights[[i_source]]=rep(1,dim(testing_input[[i_source]])[1])
    }
  }
  
  for(i_source in 1:object@num_sources){
    if(!is.null(interval_est)){
      record_interval=matrix(0,dim(testing_input[[i_source]])[1],length(interval_est[[i_source]]));
    }
    SS=floor(dim(object@post_theta)[1]/n_thinning)
    
    
    if(object@discrepancy_type[i_source]=='no-discrepancy'){
      record_cm_pred=0
      record_cm_pred_no_mean=0
      c_prop=1/4
      
      for(i_S in (1: SS)*n_thinning ){
        #print(i_S)
        
        if(i_S==floor(SS*n_thinning*c_prop)){
          cat(c_prop*100, 'percent is completed \n')
          c_prop=c_prop+1/4
        }
        
        
        theta=object@post_theta[i_S,] ##shared parameter
        sigma_2_delta=object@post_individual_par[[i_source]][i_S,1] #the first one is the sigma_2
        
        mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],math_model[[i_source]]);
        mean_cm_test_no_mean=mean_cm_test
        if(object@have_trend[i_source]){
          theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,2:(1+object@q[i_source]) ])
          
          #object@post_sample[i_S,(object@p_theta+2):(object@p_theta+1+object@q)]
          mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
        }
        
        
        record_cm_pred=record_cm_pred+mean_cm_test
        record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
        
        if(!is.null(interval_est)){
          qnorm_all=qnorm(interval_est[[i_source]]);
          for(i_int in 1:length(interval_est[[i_source]]) ){
            record_interval[,i_int]=record_interval[,i_int]+mean_cm_test+qnorm_all[i_int]*sqrt(sigma_2_delta/testing_output_weights[[i_source]])
          }
        }
        
      }
      
      record_cm_pred=record_cm_pred/SS
      record_cm_pred_no_mean=record_cm_pred_no_mean/SS
      #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_cm_pred,nrow = 64, ncol = 64,main='real')
      
      #output.list[[i_source]]=list()
      
      #output.list[[i_source]]$math_model_mean=record_cm_pred
      
      predictobj@math_model_mean[[i_source]]=record_cm_pred
      predictobj@math_model_no_mean[[i_source]]=record_cm_pred_no_mean
      if(!is.null(interval_est)){
        record_interval=record_interval/SS
        #output.list[[i_source]]$interval=record_interval
        predictobj@interval[[i_source]]=record_interval
        #ans.list[[2]]=record_interval
      }
      cat('Source',i_source, 'is completed \n')
      
      # return(output.list)
      
    }else{
      
      if(!is.null(interval_est)){
        c_star_record=rep(0,dim(testing_input[[i_source]])[1])
      }
      
      N_testing=dim(testing_input[[i_source]])[1]
      r0=as.list(1:object@p_x[i_source])
      for(i in 1:object@p_x[i_source]){
        r0[[i]]=abs(outer(object@input[[i_source]][,i],testing_input[[i_source]][,i],'-'))
      }
      
      record_cm_pred=0
      record_cm_pred_no_mean=0
      record_pred_mean=0
      
      SS=floor(dim(object@post_theta)[1]/n_thinning)
      
      #c_star=rep(0, n_testing)
      c_prop=1/4
      
      
      
      for(i_S in (1: SS)*n_thinning ){
        #print(i_S)
        
        if(i_S==floor(SS*n_thinning*c_prop)){
          cat(c_prop*100, 'percent is completed \n')
          c_prop=c_prop+1/4
        }
        
        theta=object@post_theta[i_S,]
        beta_delta=exp(object@post_individual_par[[i_source]][i_S,1:object@p_x[i_source]])
        eta_delta=exp(object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+1])
        sigma_2_delta=object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+2]
        if(object@have_trend[i_source]){
          theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
        }
        
        mean_cm=mathematical_model_eval(object@input[[i_source]],theta,object@simul_type[[i_source]],object@emulator[[i_source]],math_model[[i_source]]);
        
        output_minus_cm=object@output[[i_source]]- mean_cm
        
        if(object@have_trend[i_source]){
          output_minus_cm=output_minus_cm-object@X[[i_source]]%*%theta_m
        }
        ##remember 
        if(object@discrepancy_type[i_source]=="GaSP"){
          L=Get_R_new( beta_delta,   eta_delta,  
                       object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
                       1/object@output_weights[[i_source]])
        }else if(object@discrepancy_type[i_source]=="S-GaSP"){
          L=Get_R_z_new( beta_delta,   eta_delta,  object@tilde_lambda[i_source],
                         object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
                         1/object@output_weights[[i_source]])
        }
        
        R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm ))
        
        if(object@discrepancy_type[i_source]=="S-GaSP"){
          R_inv_y=Update_R_inv_y(R_inv_y,object@R0[[i_source]],beta_delta,object@kernel_type[i_source],object@alpha[[i_source]],
                                 object@tilde_lambda[i_source],object@num_obs[i_source])
        }
        
        r=separable_kernel(r0,beta_delta,kernel_type=object@kernel_type[i_source],object@alpha[[i_source]])
        
        rt_R_inv_y=t(r)%*%R_inv_y
        
        mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],math_model[[i_source]]);
        mean_cm_test_no_mean=mean_cm_test
        #f_M_testing=mean_cm_test
        
        if(object@have_trend[i_source]){
          mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
        }
        
        
        pred_mean=mean_cm_test+rt_R_inv_y
        
        record_cm_pred=record_cm_pred+mean_cm_test
        record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
        
        record_pred_mean=record_pred_mean+pred_mean
        
        
        
        
        if(!is.null(interval_est)){
          if(object@discrepancy_type[i_source]=="GaSP"){
            R_tilde_inv_r=(backsolve(t(L),forwardsolve(L, r )))
            
            for(i_testing in 1:N_testing){
              c_star_record[i_testing]=1-r[,i_testing]%*%R_tilde_inv_r[,i_testing]
            }
            var_gasp_f=sigma_2_delta*c_star_record
            if(interval_data[i_source]==T){
              var_gasp_f=var_gasp_f+sigma_2_delta*eta_delta/testing_output_weights[[i_source]]
            }
            
            qnorm_all=qnorm(interval_est[[i_source]]);
            for(i_int in 1:length(interval_est[[i_source]]) ){
              record_interval[,i_int]=record_interval[,i_int]+pred_mean+qnorm_all[i_int]*sqrt(var_gasp_f)
              
            }
            
          }else if(object@discrepancy_type[i_source]=="S-GaSP"){
            R=separable_kernel(object@R0[[i_source]], beta_delta, object@kernel_type[i_source],
                               object@alpha[[i_source]])
            R_tilde=R+1/object@tilde_lambda[i_source]*diag(object@num_obs[i_source]);  
            
            #system.time(Chol_Eigen(R_tilde))
            
            L_R_tilde=Chol_Eigen(R_tilde)
            
            I_minus_R_R_tilde=diag(object@num_obs[i_source])-t(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde,R )))
            
            
            R_tilde_inv_r=(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde, r )))
            
            r_z=I_minus_R_R_tilde%*%r
            
            R_z_tilde_inv_r_z=(backsolve(t(L),forwardsolve(L, r_z )))
            
            for(i_testing in 1:N_testing){
              c_star_record[i_testing]=1-t(r[,i_testing])%*%R_tilde_inv_r[,i_testing]-r_z[,i_testing]%*%R_z_tilde_inv_r_z[,i_testing]
            }
            
            t(r[,i_testing])%*%solve(R_tilde)%*%r[,i_testing] 
            
            var_gasp_f=sigma_2_delta*c_star_record
            
            if(interval_data[i_source]==T){
              var_gasp_f=var_gasp_f+sigma_2_delta*eta_delta/testing_output_weights[[i_source]]
            }
            
            qnorm_all=qnorm(interval_est[[i_source]]);
            for(i_int in 1:length(interval_est[[i_source]]) ){
              record_interval[,i_int]=record_interval[,i_int]+pred_mean+qnorm_all[i_int]*sqrt(var_gasp_f)
              
            }
            
          }
        }  
      }
      
      
      
      #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_cm_pred,nrow = 64, ncol = 64,main='real')
      
      #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_pred_mean,nrow = 64, ncol = 64,main='real')
      
      record_cm_pred=record_cm_pred/SS
      record_cm_pred_no_mean=record_cm_pred_no_mean/SS
      record_pred_mean=record_pred_mean/SS
  
      predictobj@math_model_mean[[i_source]]=record_cm_pred
      predictobj@math_model_mean_no_trend[[i_source]]=record_cm_pred_no_mean
      predictobj@mean[[i_source]]=record_pred_mean
      
      #output.list[[i_source]]=list()
      
      #output.list[[i_source]]$math_model_mean=record_cm_pred
      #output.list[[i_source]]$mean=record_pred_mean
      
      if(!is.null(interval_est)){
        record_interval=record_interval/SS
        predictobj@interval[[i_source]]=record_interval
          
        #output.list[[i_source]]$interval=record_interval
      }
      
    }
    
  }
  
  
  return(predictobj)
  
  
  #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=output.list[[1]]$interval[,2],nrow = 64, ncol = 64,main='real')
  
}