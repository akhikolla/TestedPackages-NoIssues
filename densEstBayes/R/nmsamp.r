########## R function: nmsamp ##########

# Generates a univariate normal mixture sample
# of size n for given given vectors w, mu and   
# sigma2 and a given seed.

# Last changed: 17 MAR 2004

nmsamp <- function(n,w,mu,sigma2)

{
   if (sum(w)!=1) stop("weights must sum to 1")
   if (any(sigma2<=0)) stop("variances must be positive")

   zsam <- rnorm(n)
   usam <- runif(n)
   ind <- as.vector(cut(usam,c(0,cumsum(w)),labels=F))  # Use uniform variates
   muvec <- mu[ind]                                     # to allocate component
   sigvec <- sqrt(sigma2[ind])
   
   return(sigvec*zsam + muvec)
} 
 
######## End R funtion nmsamp ########

