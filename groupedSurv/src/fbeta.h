double fbeta(double xx, void * params)
{
	double prod_l_i	= 1.0, exp_part	= 0.0, prob_dtime = 0.0; 
	exp_part = exp((*global_beta_) * global_G_[0] + xx);
	if( global_log_alpha_v_[global_Dtime_[0]-1] == -INFINITY ||
		global_log_alpha_v_[global_Dtime_[0]-1] == INFINITY)
    	prob_dtime = 0;
	else
		prob_dtime = exp(global_log_alpha_v_[global_Dtime_[0]-1]
						 *exp_part);	
	if(exp_part == INFINITY)
		prod_l_i = 0;
	else 		
		prod_l_i = likelihood_i_c(global_log_alpha_v_,
                	             global_Dtime_[0],	
								 global_Delta_[0],
								 exp_part,
								 prob_dtime);                   
	double bkb = -(1/((*global_sigma2_)*2))*xx*xx;
	double pdf_comp = exp(bkb);
	double beta_score = 0;
	// for compute score function regarding to beta			     
	if(exp_part == INFINITY || exp_part == -INFINITY)
		beta_score = 0;
    else if(exp_part <= 1e-16)
    	beta_score = 1.0;
	else	 	
		beta_score = part_likelihood_i_c(
								 global_alpha_v_,
								 global_log_alpha_v_,
                                 global_Dtime_[0],
                                 global_Delta_[0],
								 global_G_[0],
  								 exp_part,
								 prob_dtime);
	if(beta_score ==  INFINITY ||beta_score ==  -INFINITY)	
		beta_score = 0;
	double fbeta = beta_score*prod_l_i*pdf_comp;
	return fbeta;	
}


