# Wrapper function to estimate beta from family data
betaEstFam <- function(x, fam_group, fam_role, alpha, var, gtime, delta, lower, 
  upper) {
  m <- length(unique(fam_group))
  .Call("_groupedSurv_betaEst", PACKAGE = "groupedSurv", fam_group, alpha, gtime, 
    delta, x, var, lower, upper, fam_role, m)
}

