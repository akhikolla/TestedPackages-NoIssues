# function to estimate the variance of the nuisance parameter fir data with
# family struture
varEstFam <- function(x, fam_group, fam_role, alpha, gtime, delta, lower, upper, beta=0)
{
	m <- length(unique(fam_group))
  .Call("_groupedSurv_varEst", PACKAGE = "groupedSurv", fam_group, alpha, gtime, 
    delta, x, beta, lower, upper, fam_role, m)
}

