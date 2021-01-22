print.fclust <- function (x,...)
{
fclust.obj=x
if ((missing(fclust.obj)) || (!inherits(fclust.obj, "fclust")))
stop("An object of class fclust must be given")
X = if(is.null(fclust.obj$X)){fclust.obj$D}else{fclust.obj$X}
U=fclust.obj$U
k=fclust.obj$k
crit = fclust.obj$criterion
n=nrow(X)
cat("\n Fuzzy clustering object of class 'fclust' ")
cat("\n ")
cat("\n Number of objects: \n", n);
cat("\n ")
cat("\n Number of clusters: \n", k);
cat("\n ")
cat("\n Clustering index values: \n")
print(crit)
cat("\n ")
cat("\n Closest hard clustering partition: \n")
print(cl.memb(U)[,1])
cat("\n Membership degree matrix (rounded): \n")
print(round(U,2))
cat("\n Available components: \n")
print(names(fclust.obj))
cat("\n ")
invisible(fclust.obj)
}
