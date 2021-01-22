varLLFam <- function(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0)
{
	m <- length(unique(fam_group))
  .Call("_groupedSurv_ll", PACKAGE = "groupedSurv", fam_group, alpha, gtime, delta, 
    x, beta, var, fam_role, m)
}


