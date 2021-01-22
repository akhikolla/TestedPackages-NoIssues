#' Computes Randomization Metrics (explained in paper) about a design algorithm
#' 
#' @param designs	A matrix where each column is one design.
#' 
#' @return 			A list of resulting data: the probability estimates for
#' 					each pair in the design of randomness where estmates close
#' 					to ~0.5 represent random assignment, then the entropy metric
#' 					the distance metric, the maximum eigenvalue of the allocation
#' 					var-cov matrix (operator norm) and the squared Frobenius norm 
#' 					(the sum of the squared eigenvalues)
#' 
#' @author Adam Kapelner
#' @export
compute_randomization_metrics = function(designs){
	n = nrow(designs)
	r = ncol(designs)
	
	gc() #Delete at your own risk!
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("DesignMetrics.RandomizationMetrics")
	.jcall(java_obj, "V", "setNandR", as.integer(n), as.integer(r))
#	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	
	#feed in the data
	for (j in 1 : r){
		.jcall(java_obj, "V", "setDesign", as.integer(j - 1), as.integer(designs[, j])) #java indexes from 0...n-1
	}
	#get it going
	.jcall(java_obj, "V", "compute")
	
	#harvest the data and return it as a list
	p_hat_ijs = .jcall(java_obj, "[[D", "getPhats", simplify = TRUE)
	rand_entropy_metric = .jcall(java_obj, "D", "getRandEntropyMetric")
	rand_norm_se_metric = .jcall(java_obj, "D", "getRandStdErrMetric")
	
	#for the maximum eigenvalue we need to transform the allocation vector to be in {-1, 1}
	designs[designs == 0] = -1
	#now take the eigendecomposition of the variance-covariance matrix of the allocations
	e_d = eigen(var(t(designs)))	
	
	list(
		p_hat_ijs = p_hat_ijs, 
		rand_entropy_metric = rand_entropy_metric, 
		rand_norm_se_metric = rand_norm_se_metric,
		max_eigenval = max(e_d$values),
		frob_norm_sqd = sum(e_d$values^2)
	)
}