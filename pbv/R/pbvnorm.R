## File Name: pbvnorm.R
## File Version: 0.01

pbvnorm <- function(x, y, rho)
{
    res <- pbv_rcpp_pbvnorm(x=x, y=y, rho=rho)
    return(res)
}
