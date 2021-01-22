double fll(double xx, void * params)
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
    //double fll  = prod_l_i*pdf_comp/sqrt(*global_sigma2_*2*PI);
    double fll  = prod_l_i*pdf_comp;
	return fll;
}

double fsigma(double xx, void * params)
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
    double fsigma  = prod_l_i*pdf_comp*(-bkb-1)/(*global_sigma2_);
	return fsigma;
}








