eBsc.default <-
function(y, q = NULL, method = "N", parallel = FALSE,  R0 = NULL, zero_range = c(-45,1), ARpMAq = c(0,0), trace = FALSE, tol.lambda = 0.01, tol.rho = 0.01, max.iter = 50)
{
est <- EBCparallel(y, q, method, parallel, R0, zero_range, ARpMAq, trace, tol.lambda, tol.rho, max.iter)

est$q.hat      <- est$q.hat
est$lambda.hat <- est$lambda.hat
est$R.hat      <- est$R.hat
est$f.hat      <- est$f.hat
est$etq.hat    <- est$etq.hat
est$iterations <- est$iterations
est$sigma2.hat <- est$sigma2.hat
est$data       <- est$data
est$out        <- est$out

est$call <- match.call()
class(est) <- "eBsc"
est
}

