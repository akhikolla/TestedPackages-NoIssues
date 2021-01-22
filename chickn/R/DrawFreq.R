#' @title Draw frequency vectors
#' 
#' @description Function samples frequency vectors from the selected frequency distribution law.
#' @param m Number of frequency vectors.
#' @param n Length of frequency vector.
#' @param sigma Data variance, a scalar or a vector in the case of the Gaussian distribution mixture.    
#' @param alpha Variance weights. By default all are equal to 1. 
#' @param ncores Number of cores. Multicore computation should be used only when the data is a mixture of Gaussian distributions.
#' @param TypeDist Frequency distribution type. Possible values: "G" (Gaussian), 
#' "FG" (Folded Gaussian radial) or "AR" (Adapted radius). Default is "AR". 
#' @param parallel logical parameter that defines whether to perform the parallel computations. Default is FALSE.
#' @return A matrix m x n, with frequency vectors in rows.
#' @details The frequency vectors \eqn{w_1, \dots, w_m} are randomly sampled from 
#' the predefined frequency distribution. The distribution law can be either 
#' \eqn{N(0, \Sigma^{-1})} (\code{typeDist} = "G") or \eqn{p_R \cdot \varphi \cdot \Sigma^{-\frac{1}{2}}} 
#' (\code{typeDist} = c("FG", "AR")), where \eqn{\varphi} is a vector 
#' uniformly distributed on the unit sphere, \eqn{\Sigma} is a diagonal matrix with the data variance \code{sigma} on the diagonal
#' and where \eqn{p_R} is the radius density function.   
#' For "FG" the radius distribution is \eqn{N(0,1)^+} and for "AR"
#' \eqn{p_R = C \cdot (R^2 + \frac{R^4}{4})^{0.5} \cdot \exp{(-0.5 \cdot R^2)}}, where C is a normalization constant.
#' @examples 
#' W1 = DrawFreq(m = 20, n = 10, sigma = 1e-3, TypeDist = "AR")
#' W2 = DrawFreq(m = 20, n = 10, sigma = 1e-3, TypeDist = "FG")
#' W3 = DrawFreq(m = 20, n = 10, sigma = 1e-3, TypeDist = "G")
#' @seealso \code{\link{EstimSigma}}, \code{\link{GenerateFrequencies}}, \code{\link{Sketch}}
#' @note The implemented method of the frequency sampling has been proposed in \insertRef{DBLP:journals/corr/KerivenBGP16}{chickn}. 
#' @importFrom Rdpack reprompt
#' @importFrom MASS mvrnorm
#' @importFrom mvnfast rmvn
#' @importFrom zipfR Rgamma Rgamma.inv
#' @importFrom stats runif rnorm
#' @importFrom foreach %dopar%
#' @export
DrawFreq<-function(m,n, sigma, 
                   alpha = rep(1,length(sigma)),
                   TypeDist = 'AR', ncores = 1, parallel = FALSE){
  
  W = matrix(0, nrow = m, ncol = n)
  
  if (length(sigma) >1){
    cum_sum_alpha = cumsum(c(0, alpha/sum(alpha)))
    ## Generate index of Gaussian according to weghts
    IndGaussian = sapply(runif(m, min = 0, max = 1), function(x) max(which(cum_sum_alpha < x)))
  }
  switch(TypeDist,
  'G'={ # Gaussian distribution
   if(length(sigma) > 1){
     if(parallel){
       switch (Sys.info()['sysname'],
               'Linux' = {
                 cluster_Type = 'FORK'
               },
               'Windows' = {
                 cluster_Type = 'PSOCK'
               },
               "Darwin" = {
                 cluster_Type = 'FORK'
               },
               "SunOS" = {
                 cluster_Type = 'PSOCK'
               },
               stop(paste("Package is not compatible with", Sys.info()['sysname']))
       )
       cl <- parallel::makeCluster(ncores, type = cluster_Type)
       doParallel::registerDoParallel(cl)
       doRNG::registerDoRNG()
       on.exit(parallel::stopCluster(cl))
       W <- foreach::foreach(j = 1:m, .combine = rbind) %dopar%{
         mvrnorm(n=1, mu = rep(0, n), Sigma = solve(sigma[IndGaussian[j]]*diag(n)))
       }
     }else{
       for(j in 1:m){
         W[j,]=mvrnorm(n=1, mu = rep(0, n), Sigma = solve(sigma[IndGaussian[j]]*diag(n)))
       }
     }
   }else{
     W = mvrnorm(n=m, mu = rep(0, n), Sigma = sigma^(-1)*diag(n))
   }
  },
  'FG'={
   R = abs(rnorm(m)) # Folded Gaussian distribution
  },
  'AR' = {
   ## R ~ p_R = C*(R^2 + R^4/4)^1/2*exp(-1/2*R^2)
   ## C = sqrt(2)/exp(2)*gamma(1.5)*gammainc(2,1.5)
   Funx_inv = function(x) 2*sqrt(x/2 - 1)
   Funy = function(y) y*(1 - Rgamma(a = 1.5, x = 2)) + Rgamma(a = 1.5, x = 2)
   R = sapply(runif(m, min = 0, max = 1), function(x) Funx_inv(Rgamma.inv(a = 1.5, y = Funy(x))))
  },
  stop('Distribution type is unknown'))
  if(TypeDist != 'G'){
    #Draw directions on unit shere from uniform distribution
    phi = matrix(0, nrow = m, ncol = n)
    phi = rmvn(n=m, mu = rep(0, n), sigma = diag(n), ncores = ncores)
    phi = phi/sqrt(rowSums(phi*phi))
    if(length(sigma) == 1){
      InvSigma = rep(sigma^-0.5, n)
      W = R*matrix(InvSigma, nrow = m, ncol = n, byrow = TRUE)*phi
    }else{
      if(parallel){
        switch (Sys.info()['sysname'],
                'Linux' = {
                  cluster_Type = 'FORK'
                },
                'Windows' = {
                  cluster_Type = 'PSOCK'
                },
                "Darwin" = {
                  cluster_Type = 'FORK'
                },
                "SunOS" = {
                  cluster_Type = 'PSOCK'
                },
                stop(paste("Package is not compatible with", Sys.info()['sysname']))
        )
        cl <- parallel::makeCluster(ncores, type = cluster_Type)
        doParallel::registerDoParallel(cl)
        doRNG::registerDoRNG()
        on.exit(parallel::stopCluster(cl))
        W <- foreach::foreach (j = 1:m, .combine = rbind) %dopar%{
          R[j]*phi[j,]%*%solve(sigma[IndGaussian[j]]*diag(n))^0.5
        }
      }else{
        for(j in 1:m){
          W[j,] = R[j]*phi[j,]%*%solve(sigma[IndGaussian[j]]*diag(n))^0.5
        }
      }
    }
  }
  return(W)
}
