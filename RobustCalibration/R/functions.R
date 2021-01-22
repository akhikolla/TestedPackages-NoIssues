

mathematical_model_eval<-function(input,theta,simul_type, emulator,math_model){
  #VectorXd cm_obs_cur;
  p_theta=length(theta);
  p_x=dim(input)[2];
  n_testing=dim(input)[1];
  
  if(simul_type==0){ ###emulator
    testing_input=cbind(input, t(matrix(theta,p_theta,n_testing)));
    cm_obs_cur=Sample(emulator,testing_input);
  }else if(simul_type==1){
    testing_input=cbind(input, t(matrix(theta,p_theta,n_testing)));
    cm_obs_cur=math_model(input, theta);
  }else if( (simul_type==2) | (simul_type==3) ){
    cm_obs_cur= Mogihammer(input,theta,simul_type);
  }
  return(cm_obs_cur);
}


##tell whether it has discrepancy function or no
post_sample <- function(input, output, R0_list, kernel_type, p_theta, output_weights,
                        par_cur, tilde_lambda,prior_par,theta_range,S, X, have_trend, alpha,sd_proposal,
                        discrepancy_type, simul_type,emulator,math_model){ 
  #method_type=Input_list[[17]];
  if(discrepancy_type=='no-discrepancy'){
    ans=post_sample_no_discrepancy(input, output, R0_list,  p_theta, output_weights,
                                   par_cur, theta_range,S, X, have_trend, alpha,sd_proposal,
                                   discrepancy_type, simul_type,emulator,math_model);
  }else{
    ans=post_sample_with_discrepancy(input, output, R0_list, kernel_type, p_theta, output_weights,
                                     par_cur, tilde_lambda,prior_par,theta_range,S, X, have_trend, alpha,sd_proposal,
                                     discrepancy_type, simul_type,emulator,math_model);
  }
  return(ans);
  
}

##sample the posterior with discrepancy
post_sample_with_discrepancy<-function(input, output, R0_list, kernel_type, p_theta, output_weights,
                             par_cur, tilde_lambda,prior_par,theta_range,S, X, have_trend, alpha,sd_proposal,
                             discrepancy_type,simul_type,emulator,math_model){
  

  mylist=as.list(1:4);
  
  #Initialization
  #MatrixXd input=Input_list[0];
  num_obs=dim(input)[1];
  p_x=dim(input)[2];
  
  #MatrixXd output=Input_list[1];
  #List R0_list=Input_list[2];
  #String   kernel_type=Input_list[3];
  #int      p_theta=Input_list[4];
  #VectorXd area_size=Input_list[6];
  #VectorXd par_cur=Input_list[7]; //initial values
  #double   tilde_lambda=0;  //here if it is GaSP it means tilde_lambda=0
  #int method_type=Input_list[16];  //if it is GaSP, it is 1. if it is S-GaSP, it is 2.
  #cm_eval=simul_type; #this is the computer model
  #Input_list[17]; //this is the computer model
  
  #if(method_type==2){
  #  tilde_lambda=Input_list[8];
  #}
  CL_a_b=prior_par;
  
  CL=CL_a_b[1:p_x];
  a=CL_a_b[p_x+1];
  b=CL_a_b[p_x+2];
  
  #MatrixXd theta_range=Input_list[10];
  #int      S=Input_list[11];
  
  #VectorXd X=Input_list[12];
  #have_mean=have_trend;
  p_theta_m=0;
  
  if(have_trend){
    p_theta_m=dim(X)[2];
    #.cols();
  }
  
  #VectorXd alpha=  Input_list[14];
  #VectorXd sd_all=  Input_list[15];
   sd_theta=sd_proposal[1:p_theta];
   sd_log_beta=sd_proposal[(p_theta+1):(p_theta+p_x)];
   sd_log_eta=sd_proposal[p_theta+p_x+1];
  
  
  inv_output_weights=1/output_weights;
  

  record_par=matrix(0,S,p_theta+p_x+2+p_theta_m);
  record_post=rep(0,S)
  #MatrixXd::Zero(S,p_theta+p_x+2+p_theta_m);
  #VectorXd record_post=VectorXd::Zero(S);
  
  param=par_cur;
  
  
  #MatrixXd L;
  if(discrepancy_type=='GaSP'){
    L=Get_R_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),R0_list,kernel_type,alpha, inv_output_weights);
  }else{
    L=Get_R_z_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),tilde_lambda,R0_list,kernel_type,alpha, inv_output_weights );
  }
  
  #MatrixXd L_propose;
  
  accept_S_theta=0;
  accept_S_beta=0;
  count_dec_record=rep(0,S); # this is to count how many sample points are outside the boundary and get rejected
  
  post_cur=0;
  
  theta_cur=rep(0,p_theta);
  theta_sample=rep(0,p_theta);
  
  #bool decision_0;
  #bool decision;
  param_propose=par_cur;
  xi_sample=rep(0,p_x);
  
  log_eta_sample=0;
  #post_propose;
  r_ratio=0;


  cm_obs_cur=mathematical_model_eval(input,param[1:p_theta],simul_type,emulator,math_model);
  
  
  c_prop=1/4
  
  #start of the sampling
  
  for (i_S in 1:S){
    if(i_S==floor(S*c_prop)){
      #cat(post_cur,'\n')
      cat(c_prop*S, ' of ', S, ' posterior samples are drawn \n')
      c_prop=c_prop+1/4
    }
    #cat(post_cur,'\n')
    #cat(par_cur,'\n')
    #cat(accept_S_theta,'\n')
    
    par_cur[(p_theta+p_x+2):(p_theta+p_x+2+p_theta_m)]=Sample_sigma_2_theta_m(par_cur,L,output,p_theta,p_x,X,have_trend,cm_obs_cur);
    
   # par_cur.segment(p_theta+p_x+1,1+p_theta_m)=Sample_sigma_2_theta_m(par_cur,L,output,p_theta,p_x,X,have_trend,cm_obs_cur);
    
    post_cur=Log_marginal_post(par_cur,L,output,p_theta,p_x,X,have_trend,CL,a,b,cm_obs_cur);
    
    #//sample theta
    theta_cur=par_cur[1:p_theta];
    decision_0=F;
    
    
    # for(int i_theta=0; i_theta<p_theta; ++i_theta){
    #   theta_sample(i_theta)=theta_cur(i_theta)+sd_theta(i_theta)*(theta_range(i_theta,1)-theta_range(i_theta,0))*distribution_stan_norm(generator);
    #   if((theta_sample(i_theta)>theta_range(i_theta,1))|| (theta_sample(i_theta)<theta_range(i_theta,0)) ){
    #     decision_0=true;
    #     count_dec_record(i_S)=1;
    #     break;
    #   }
    #   
    # }
    
    
    for(i_theta in 1:p_theta){
      theta_sample[i_theta]=rnorm(1,mean=theta_cur[i_theta],sd= (sd_theta[i_theta]*(theta_range[i_theta,2]- theta_range[i_theta,1]) ) )  ##here maybe truncated
       #cat(theta_sample[i_theta],'\n')
       #cat(theta_range[i_theta,2],'\n')
       #cat(decision_0,'\n')
        if((theta_sample[i_theta]>theta_range[i_theta,2])| (theta_sample[i_theta]<theta_range[i_theta,1]) ){
          decision_0=T  ##reject directly
        count_dec_record[i_S]=1
        break;
      }
    }


    #//if decision_0==true, then stay at original place
    #//ow see whether we accept the sample
    if(decision_0==F){
      param_propose=par_cur;
      param_propose[1:p_theta]=theta_sample;
      

      cm_obs_propose=mathematical_model_eval(input,theta_sample,simul_type,emulator,math_model);
      
      post_propose=Log_marginal_post(param_propose,L,output,p_theta,p_x,X,have_trend,CL,a,b, cm_obs_propose);
      r_ratio=exp(post_propose-post_cur);
      decision=Accept_proposal(r_ratio);
      if(decision){

        par_cur=param_propose;
        post_cur=post_propose;
        cm_obs_cur=cm_obs_propose;
        accept_S_theta=accept_S_theta+1;
        
        #cat(par_cur,'\n')
        
      }
      
    }
    
    #//sample xi and log_eta 
    for(i_x in 1: p_x){
      xi_sample[i_x]=par_cur[p_theta+i_x]+sd_log_beta[i_x]*rnorm(1);
      #distribution_stan_norm(generator);
    }
    log_eta_sample=par_cur[p_theta+p_x+1]+sd_log_eta*rnorm(1); 
    #distribution_stan_norm(generator);
    
    param_propose=par_cur;
    
    param_propose[(p_theta+1):(p_theta+p_x)]=xi_sample;
    param_propose[p_theta+p_x+1]=log_eta_sample;
    
    param=param_propose;
    
    
    if(discrepancy_type=='GaSP'){
      L_propose=Get_R_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),R0_list,kernel_type,alpha, inv_output_weights);
    }else{
      L_propose=Get_R_z_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),tilde_lambda,R0_list,kernel_type,alpha, inv_output_weights );
    }
    
    post_propose=Log_marginal_post(param,L_propose,output,p_theta,p_x,X,have_trend,CL,a,b,cm_obs_cur);
    
    r_ratio=exp(post_propose-post_cur);
    
    decision=Accept_proposal(r_ratio);
    
    if(decision){
      par_cur=param_propose;
      post_cur=post_propose;
      L=L_propose;
      accept_S_beta=accept_S_beta+1;
    }
    
    record_par[i_S,]=par_cur;
    record_post[i_S]=post_cur;
    
    #record_par.block(i_S,0,1,p_theta+p_x+2+p_theta_m)=par_cur.transpose();
    #record_post(i_S)=post_cur;
  }
  

  mylist[[1]]=record_par;
  mylist[[2]]=record_post;
  
  mylist[[3]]=c(accept_S_theta,accept_S_beta);
  
  mylist[[4]]=count_dec_record;
  
  cat('Done with posterior sampling \n')
  
  cat(accept_S_theta, ' of ', S, ' proposed posterior samples of calibration parameters are accepted \n')
  cat(accept_S_beta, ' of ', S, ' proposed posterior samples of range and nugget parameters are accepted \n')
  
  return(mylist);
  
  
}






post_sample_no_discrepancy<-function(input, output, R0_list,  p_theta, output_weights,
                                par_cur, theta_range,S, X, have_trend, alpha,sd_proposal,
                                discrepancy_type, simul_type,emulator,math_model){
  
  
  mylist=list();
  
  #Initialization
  num_obs=dim(input)[1];
  p_x=dim(input)[2];

  p_theta_m=0;
  
  if(have_trend){
    p_theta_m=dim(X)[2];
  }
  
  #List mylist;
  sd_theta=sd_proposal[1:p_theta];
  
  # simul_type=  Input_list[17];
  

  
  inv_output_weights=1/output_weights;

  
  record_par=matrix(0,S,p_theta+1+p_theta_m);
  record_post=rep(0,S)

  accept_S_theta=0;
  count_dec_record=rep(0,S); # this is to count how many sample points are outside the boundary and get rejected
  
  
  #MatrixXd record_par=MatrixXd::Zero(S,p_theta+1+p_theta_m);
  #VectorXd record_post=VectorXd::Zero(S);
  
  param=par_cur;
  

  post_cur=0;
  
  
  theta_cur=rep(0,p_theta);
  theta_sample=rep(0,p_theta);
  
  #bool decision_0;
  #bool decision;
  param_propose=par_cur;
  r_ratio=0;

  
  cm_obs_cur=mathematical_model_eval(input,param[1:p_theta],simul_type,emulator,math_model);
  
  
  #cm_obs_cur=mathematical_model_eval(param.head(p_theta),input,simul_type,emulator_par,simul_input,simul_output);
  
  c_prop=1/4
  for (i_S in 1:S){
    if(i_S==floor(S*c_prop)){
      #cat(post_cur,'\n')
      cat(c_prop*S, ' of ', S, ' posterior samples are drawn \n')
      c_prop=c_prop+1/4
    }
    #cat(post_cur,'\n')
    
    #cat(par_cur,'\n')
    par_cur[(p_theta+1):(p_theta+1+p_theta_m)]=Sample_sigma_2_theta_m_no_discrepancy(par_cur,output,p_theta,X,have_trend, inv_output_weights,cm_obs_cur);
    
    # par_cur.segment(p_theta+p_x+1,1+p_theta_m)=Sample_sigma_2_theta_m(par_cur,L,output,p_theta,p_x,X,have_trend,cm_obs_cur);
    
    post_cur=Log_marginal_post_no_discrepancy(par_cur,output,p_theta,X,have_trend, inv_output_weights,cm_obs_cur);
    

    #par_cur.segment(p_theta,1+p_theta_m)=Sample_sigma_2_theta_m_no_discrepancy(par_cur,output,p_theta,X,have_mean, inv_output_weights,cm_obs_cur);
    
    #post_cur=Log_marginal_post_no_discrepancy(par_cur,output,p_theta,X,have_mean, inv_output_weights,cm_obs_cur);
    
    #sample theta
    theta_cur=par_cur[1:p_theta];
    decision_0=F;
    
    for(i_theta in 1:p_theta){
      theta_sample[i_theta]=rnorm(1,mean=theta_cur[i_theta],sd=sd_theta[i_theta]*(theta_range[i_theta,2]- theta_range[i_theta,1]) )  ##here maybe truncated
      if((theta_sample[i_theta]>theta_range[i_theta,2])| (theta_sample[i_theta]<theta_range[i_theta,1]) ){
        decision_0=T  ##reject directly
        count_dec_record[i_S]=1
        break;
      }
    }
    
    #//if decision_0==true, then stay at original place
    #//ow see whether we accept the sample

    if(decision_0==F){
      param_propose=par_cur;
      param_propose[1:p_theta]=theta_sample;
      cm_obs_propose=mathematical_model_eval(input,theta_sample,simul_type,emulator,math_model);
      
      post_propose=Log_marginal_post_no_discrepancy(param_propose,output,p_theta,X,have_trend, inv_output_weights,cm_obs_propose);

      #  Log_marginal_post(param_propose,L,output,p_theta,p_x,X,have_trend,CL,a,b, cm_obs_propose);
      r_ratio=exp(post_propose-post_cur);
      decision=Accept_proposal(r_ratio);
      if(decision){
        par_cur=param_propose;
        post_cur=post_propose;
        cm_obs_cur=cm_obs_propose;
        accept_S_theta=accept_S_theta+1;
      }
      
    }
    

    
    record_par[i_S,]=par_cur;
    record_post[i_S]=post_cur;
    
    #record_par.block(i_S,0,1,p_theta+p_x+2+p_theta_m)=par_cur.transpose();
    #record_post(i_S)=post_cur;
  }
  
  cat('Done with posterior sampling \n')
  mylist[[1]]=record_par;
  mylist[[2]]=record_post;
  
  mylist[[3]]=c(accept_S_theta);
  
  mylist[[4]]=count_dec_record;
  
  cat(accept_S_theta, ' of ', S, ' proposed posterior samples samples of calibration parameters are accepted \n')
  
    
  return(mylist);
  
  
}




post_sample_MS <- function(model,par_cur_theta, par_cur_individual, emulator,math_model_MS){ 
  
  mylist=list();

  ##################
  #3177.183 2593.284
  #[1] -6.474958e+00 -6.677884e+00 -1.070080e+01  1.854997e-05  8.482068e-03
  #[1] 3.267227e+01 7.877541e+02 1.536619e+03 1.464163e-02 3.003662e-01
  #[1] 3186.917 2595.137
  
  record_par_individual=as.list(rep(0,model@num_sources))
  for(i in 1: model@num_sources){
    if(model@discrepancy_type[i]=='GaSP' | model@discrepancy_type[i]=='S-GaSP'){
      record_par_individual[[i]]=  matrix(0,model@S,model@p_x[i]+2+model@q[i]);
    }else{
      record_par_individual[[i]]=  matrix(0,model@S,1+model@q[i]);
      
    }
  }
  
  record_post=matrix(0,model@S,model@num_sources)
  record_theta=matrix(0,model@S,model@p_theta)
  
  
  #MatrixXd L;
  L=as.list(rep(0,model@num_sources))
  L_propose=L
  
  for(i in 1: model@num_sources){
      if(model@discrepancy_type[i]=='GaSP'){
        L[[i]]=Get_R_new(exp(par_cur_individual[[i]][1:model@p_x[i]]),exp(par_cur_individual[[i]][model@p_x[i]+1]),
                         model@R0[[i]],model@kernel_type[i],model@alpha[[i]], 1/model@output_weights[[i]]);
      }else if(model@discrepancy_type[i]=='S-GaSP'){

        L[[i]]=Get_R_z_new(exp(par_cur_individual[[i]][1:model@p_x[i]]),exp(par_cur_individual[[i]][model@p_x[i]+1]),
                           model@tilde_lambda[i],model@R0[[i]],model@kernel_type[i],
                           model@alpha[[i]], 1/model@output_weights[[i]] );
      }
  }
  #L_propose=L
  #MatrixXd L_propose;
  
  accept_S_theta=0;
  accept_N_cov_par=rep(0,model@num_sources);
  accept_S_beta=rep(0,model@num_sources);
  
  count_dec_record=rep(0,model@S); # this is to count how many sample points are outside the boundary and get rejected for calibration
  
  #post_cur=0;
  
  theta_cur=par_cur_theta;
  theta_sample=rep(0,model@p_theta);
  
  individual_par_cur=par_cur_individual
  
  ###########################
  #param_propose=par_cur;
  #xi_sample=rep(0,p_x);
  #log_eta_sample=0;
  #post_propose;
  r_ratio=0;
  
  cm_obs_cur=as.list(rep(0,model@num_sources))
  cm_obs_propose=cm_obs_cur
  
  for(i in 1:model@num_sources){
    cm_obs_cur[[i]]=mathematical_model_eval(model@input[[i]],theta_cur[model@index_theta[[i]]],model@simul_type[i],emulator[[i]],math_model_MS[[i]]);
  }
  
  c_prop=1/4
  
  
  post_cur=rep(0,model@num_sources)
  post_propose=rep(0,model@num_sources)
  
  for (i_S in 1:model@S){
    if(i_S==floor(model@S*c_prop)){
      #cat(post_cur,'\n')
      cat(c_prop*model@S, ' of ', model@S, ' posterior samples are drawn \n')
      c_prop=c_prop+1/4
    }

    #print(i_S)
    #print(accept_S_theta)
    #print(individual_par_cur[[1]])
    #print(theta_cur)
    #print(post_cur)
    #print(individual_par_cur[[2]])
    
    #print(post_cur)
    
    
    for(i in 1:model@num_sources){
      if(model@discrepancy_type[i]=='no-discrepancy'){
        sigma_2_theta_m_sample=Sample_sigma_2_theta_m_no_discrepancy(c(theta_cur[model@index_theta[[i]]], 
                                                        individual_par_cur[[i]]),
                                                      model@output[[i]],length(model@index_theta[[i]]),
                                                      model@X[[i]],
                                                      model@have_trend[i],1/model@output_weights[[i]],cm_obs_cur[[i]] );
        
        
        individual_par_cur[[i]]=sigma_2_theta_m_sample
        
        post_cur[i]=Log_marginal_post_no_discrepancy(c(theta_cur[model@index_theta[[i]]], 
                                                       individual_par_cur[[i]]),model@output[[i]],length(model@index_theta[[i]]),
                                                     model@X[[i]],
                                                     model@have_trend[i],1/model@output_weights[[i]],cm_obs_cur[[i]] );
        record_par_individual[[i]][i_S,]=sigma_2_theta_m_sample  ##first par is the sigma_2
      }else{
        sigma_2_theta_m_sample=Sample_sigma_2_theta_m(c(theta_cur[model@index_theta[[i]]], 
                                                        individual_par_cur[[i]]),L[[i]],
                                                      model@output[[i]],
                                                      length(model@index_theta[[i]]),
                                                      model@p_x[i],model@X[[i]],
                                                      model@have_trend[i],cm_obs_cur[[i]] );
        
        
        individual_par_cur[[i]][ (model@p_x[i]+2):(model@p_x[i]+2+model@q[i]) ]=sigma_2_theta_m_sample
        
        post_cur[i]=Log_marginal_post(c(theta_cur[model@index_theta[[i]]], 
                                        individual_par_cur[[i]]),L[[i]],
                                      model@output[[i]],
                                      length(model@index_theta[[i]]),
                                      model@p_x[i],model@X[[i]],
                                      model@have_trend[i],
                                      model@prior_par[[i]][1:model@p_x[i]],
                                      model@prior_par[[i]][model@p_x[i]+1],
                                      model@prior_par[[i]][model@p_x[i]+2],
                                      cm_obs_cur[[i]] );
        
        
        
        
        #//sample xi and log_eta 
        xi_sample=rep(0,model@p_x[i])
        
        for(i_x in 1: model@p_x[i]){
          xi_sample[i_x]=individual_par_cur[[i]][i_x]+model@sd_proposal_cov_par[[i]][i_x]*rnorm(1);
          #distribution_stan_norm(generator);
        }
        log_eta_sample=individual_par_cur[[i]][model@p_x[i]+1]+
                       model@sd_proposal_cov_par[[i]][model@p_x[i]+1]*rnorm(1); 
        
        individual_par_propose=  individual_par_cur[[i]]
        individual_par_propose[1:model@p_x[i]]=xi_sample
        individual_par_propose[model@p_x[i]+1]=log_eta_sample
        
        
        
        if(model@discrepancy_type[i]=='GaSP'){
          L_propose[[i]]=Get_R_new(exp(individual_par_propose[1:model@p_x[i]]),
                                   exp(individual_par_propose[model@p_x[i]+1]),
                           model@R0[[i]],model@kernel_type[i],model@alpha[[i]], 1/model@output_weights[[i]]);
        }else if(model@discrepancy_type[i]=='S-GaSP'){
          L_propose[[i]]=Get_R_z_new(exp(individual_par_propose[1:model@p_x[i]]),
                                     exp(individual_par_propose[model@p_x[i]+1]),
                             model@tilde_lambda[i],model@R0[[i]],model@kernel_type[i],
                             model@alpha[[i]], 1/model@output_weights[[i]] );
        }
        
        
        post_propose[i]=Log_marginal_post(c(theta_cur[model@index_theta[[i]]], 
                                            individual_par_propose),L_propose[[i]],
                                      model@output[[i]],
                                      length(model@index_theta[[i]]),
                                      model@p_x[i],model@X[[i]],
                                      model@have_trend[i],
                                      model@prior_par[[i]][1:model@p_x[i]],
                                      model@prior_par[[i]][model@p_x[i]+1],
                                      model@prior_par[[i]][model@p_x[i]+2],
                                      cm_obs_cur[[i]] );
        

        r_ratio=exp( post_propose[i]-post_cur[i]);
        
        decision=Accept_proposal(r_ratio);
        
        if(decision){
          individual_par_cur[[i]]=individual_par_propose;
          post_cur[i]=post_propose[i];
          L[[i]]=L_propose[[i]];
          accept_S_beta[i]=accept_S_beta[i]+1;
        }
        
        record_par_individual[[i]][i_S,]=individual_par_cur[[i]]
      }
      
    }
      
    ############################################################
    #//sample theta
    #theta_cur=par_cur[1:p_theta];
    decision_0=F;
    
    

    
    for(i_theta in 1:model@p_theta){
      theta_sample[i_theta]=rnorm(1,mean=theta_cur[i_theta],sd= (model@sd_proposal_theta[i_theta]*(model@theta_range[i_theta,2]- model@theta_range[i_theta,1]) ) )  ##here maybe truncated
      #cat(theta_sample[i_theta],'\n')
      #cat(theta_range[i_theta,2],'\n')
      #cat(decision_0,'\n')
      if((theta_sample[i_theta]>model@theta_range[i_theta,2])| (theta_sample[i_theta]<model@theta_range[i_theta,1]) ){
        decision_0=T  ##reject directly
        count_dec_record[i_S]=1
        break;
      }
    }
    
    
    #//if decision_0==true, then stay at original place
    #//ow see whether we accept the sample
    if(decision_0==F){
      #theta_propose=theta_sample;
      
      for(i in 1:model@num_sources){
        cm_obs_propose[[i]]=mathematical_model_eval(model@input[[i]],theta_sample[model@index_theta[[i]]],model@simul_type[i],emulator[[i]],math_model_MS[[i]]);

        if(model@discrepancy_type[i]=='no-discrepancy'){
          post_propose[i]=Log_marginal_post_no_discrepancy(c(theta_sample[model@index_theta[[i]]], 
                                                         individual_par_cur[[i]]),model@output[[i]],length(model@index_theta[[i]]),
                                                         model@X[[i]],
                                                         model@have_trend[i],1/model@output_weights[[i]],cm_obs_propose[[i]] );
        }else{
          post_propose[i]=Log_marginal_post(c(theta_sample[model@index_theta[[i]]], 
                                           individual_par_cur[[i]]),L[[i]],
                                         model@output[[i]],
                                         length(model@index_theta[[i]]),
                                         model@p_x[i],model@X[[i]],
                                         model@have_trend[i],
                                         model@prior_par[[i]][1:model@p_x[i]],
                                         model@prior_par[[i]][model@p_x[i]+1],
                                         model@prior_par[[i]][model@p_x[i]+2],
                                         cm_obs_propose[[i]] );
          
        }
        

      }
      
      
      r_ratio=exp(sum(post_propose)-sum(post_cur));
      decision=Accept_proposal(r_ratio);
      if(decision){
        
        theta_cur=theta_sample;
        post_cur=post_propose;
        cm_obs_cur=cm_obs_propose;
        accept_S_theta=accept_S_theta+1;
        
        #cat(par_cur,'\n')
        
      }
      
    }
    
    record_post[i_S,]=post_cur
    
    record_theta[i_S,]=theta_cur
    
    
    
    
  }
    
  mylist$record_post=record_post
    
  mylist$record_theta=record_theta
  
  mylist$individual_par=record_par_individual
  
  mylist$accept_S_theta=accept_S_theta
  mylist$accept_S_beta=accept_S_beta
  mylist$count_dec_record=count_dec_record
  
  
  cat('Done with posterior sampling \n')
  cat(accept_S_theta, ' of ', model@S, ' proposed posterior samples of calibration parameters are accepted \n')
  
  for(i in 1:model@num_sources){
    if(model@discrepancy_type[i]!='no-discrepancy'){
       cat('For source', i,':', accept_S_beta[i], ' of ', model@S, 
           ' proposed posterior samples of range and nugget parameters are accepted \n')
    }
  }
  return(mylist)
  
}
  
  
  
  
