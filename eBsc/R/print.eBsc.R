print.eBsc <-
function(x,...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("\n Empirical Bayesian Penalized Smopthing Splines with Correlated Errors:\n")
cat("\n")
cat("all parameters computed  iteratively\n")
cat("\n")
cat("\n Summary statistics of errors:\n")
print(summary(x$res))

}
