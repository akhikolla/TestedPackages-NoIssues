double falpha(double xx, void * params)
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

    int i = *global_alpha_index_;
	// *global_alphaIndex_ is the index for alpha, 
	// starting from 0, since this is defined in c++.
	double alpha_score 		= 0.0;
	if(exp_part == INFINITY)
		alpha_score = 0;
	else{
		// case 1: i = dtime - 1 
			// sub case 1: delta == 0
				//	return 0
			// sub case 2: delta == 1
				// sub sub case alpha_j = 0
					// return 0
				// sub sub case alpha_j ~=0
		if(global_Dtime_[0] == i+1 && global_Delta_[0] != 0)
		{	
			if(exp_part == 0 || exp_part == INFINITY || global_alpha_v_[i]==0)	
				alpha_score = 0;
			else
				alpha_score = -exp_part*prob_dtime/(1-prob_dtime)/global_alpha_v_[i];
		}
		// case 2: i < dtime - 1
		if(global_Dtime_[0] > i+1) 
			alpha_score = exp_part/global_alpha_v_[i]; 
		// case 3: i > dtime - 1		
		if(global_Dtime_[0] < i+1)
			alpha_score = 0;
	}
	double falpha = alpha_score*prod_l_i*pdf_comp;
	return falpha;	
}


