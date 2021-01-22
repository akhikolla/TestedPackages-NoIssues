## File Name: dbvnorm.R
## File Version: 0.03

dbvnorm <- function(x, y, rho, log=FALSE)
{
    res <- pbv_rcpp_dbvnorm(x=x, y=y, rho=rho, use_log=log)
    return(res)
}
