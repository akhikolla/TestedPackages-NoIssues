# #SimCoxIntervalCens create data set with right-censored event times and interval-censored 
# # covariates
# # Time scale is such that most events occuer in the first 5 years
# # n.sample: Scalar, sample size.
# # lambda: Paramter of the baseline hazard 
# # beta0: log hazard ratio of X, the exposure/treatment
# # mu - mean of censoring time
SimCoxIntervalCensSingle <- function(n.sample, lambda = 0.25, alpha= NULL, beta0 = log(0.75), 
                               mu = 1, n.points = 5, weib.shape = 1, weib.scale = 1)
{
  ## Simulate changepoint time
  v <- stats::rweibull(n.sample, shape = weib.shape, scale = weib.scale)
 # v <- rnorm(n.sample, 3 , 0.5)
#v <- runif(n.sample,0,10)
   #####If alpha is null, Create a piecewise exponential where the intervals are [0,V) and [V,\infty] 
  u <- stats::runif(n.sample) # Random uniform variable
  event.tm <- vector(length = n.sample)
  if (is.null(alpha))
  {
    tm.cond <- -log(u) <  lambda*v # the probability of the condition is the probability of event in the second function piece
    event.tm[tm.cond==1] <- -log(u[tm.cond==1])/lambda
    event.tm[tm.cond==0] <- -(log(u[tm.cond==0]) + lambda*v[tm.cond==0]*(1-exp(beta0)))/(lambda*exp(beta0))
    
    #event.tm <- ifelse(tm.cond, -log(u)/lambda,  -(log(u) + lambda*v*(1-exp(beta0)))/(lambda*exp(beta0)))
  } else
    ##### Create a piecewise Gompertz distribution where the intervals are [0,V) and [V,\infty] (See Austin, 2011)
    
  {
    tm.cond <- -log(u)  < (lambda/alpha)*(exp(alpha*v)-1)
    event.tm[tm.cond==1] <- (1/alpha)*log(1 + (alpha/lambda)*(-log(u[tm.cond==1])) )
    event.tm[tm.cond==0] <- (1/alpha)*log((alpha/(lambda*exp(beta0))*(-log(u[tm.cond==0])))- 
                                            (exp(alpha*v[tm.cond==0])-1-exp(beta0+alpha*v[tm.cond==0]))/(exp(beta0))) 
    #event.tm <- ifelse(tm.cond, (1/alpha)*log(1 + (alpha/lambda)*(-log(u)) ),
    #                   (1/alpha)*log((alpha/(lambda*exp(beta0))*(-log(u)))- (exp(alpha*v)-1-exp(beta0+alpha*v))/(exp(beta0)))) 
  }
    ##### Exponential Censoring + maximu followup time of 5 years
  cens <- pmin(stats::rexp(n.sample,mu), 5)
 obs.tm <- pmin(event.tm,cens)
 delta <- event.tm < cens
# obs.tm <- event.tm
#  delta <- rep(1, n.sample)
   ### Interval censoring of V ### Create questionairre times
  w <- matrix(nrow = n.sample, ncol = n.points)
  points <- seq(0,5, length.out = n.points+1)
#   intervals <- seq(0,max(event.tm), length.out = n.intervals+1)
#  intervals <- seq(0,max(event.tm), by = 5)
     for (j in 1:n.points)
  {
    w[,j] <- stats::runif(n = n.sample, min = points[j], max = points[j+1])
    #   w[,j] <- intervals[j+1]
  }
  
  # w.res is the vector of zeroes and ones (X_i(t)) where timings are the values in w
  w.res <- apply(w,2,function(x) x>v)
  w <- apply(w, 2, function(x) ifelse(x < obs.tm, x , Inf)) # Correct questionairre times after event
  w.res[w==Inf] <- Inf # Correct questionairre times after event
  
  list.back <- list(event.tm = event.tm,obs.tm = obs.tm, delta = delta, 
                    w = w, w.res = w.res, v = v)
  
  return(list.back)
}


SimCoxIntervalCensSingleNormal <- function(n.sample, lambda = 0.25, alpha= NULL, beta0 = log(0.75), 
                                     mu = 1, n.points = 5, nor.mean = 5, nor.sd = 2)
{
  ## Simulate changepoint time
  #v <- rweibull(n.sample, shape = weib.shape, scale = weib.scale)
  v <- stats::rnorm(n.sample, nor.mean , nor.sd)
  v <- ifelse(v < 0, 0, v)
  v <- ifelse(v > 10, 10, v)
  #v <- runif(n.sample,0,10)
  #####If alpha is null, Create a piecewise exponential where the intervals are [0,V) and [V,\infty] 
  u <- stats::runif(n.sample) # Random uniform variable
  event.tm <- vector(length = n.sample)
  if (is.null(alpha))
  {
    tm.cond <- -log(u) <  lambda*v # the probability of the condition is the probability of event in the second function piece
    event.tm[tm.cond==1] <- -log(u[tm.cond==1])/lambda
    event.tm[tm.cond==0] <- -(log(u[tm.cond==0]) + lambda*v[tm.cond==0]*(1-exp(beta0)))/(lambda*exp(beta0))
    
    #event.tm <- ifelse(tm.cond, -log(u)/lambda,  -(log(u) + lambda*v*(1-exp(beta0)))/(lambda*exp(beta0)))
  } else
    ##### Create a piecewise Gompertz distribution where the intervals are [0,V) and [V,\infty] (See Austin, 2011)
    
  {
    tm.cond <- -log(u)  < (lambda/alpha)*(exp(alpha*v)-1)
    event.tm[tm.cond==1] <- (1/alpha)*log(1 + (alpha/lambda)*(-log(u[tm.cond==1])) )
    event.tm[tm.cond==0] <- (1/alpha)*log((alpha/(lambda*exp(beta0))*(-log(u[tm.cond==0])))- 
                                            (exp(alpha*v[tm.cond==0])-1-exp(beta0+alpha*v[tm.cond==0]))/(exp(beta0))) 
    #event.tm <- ifelse(tm.cond, (1/alpha)*log(1 + (alpha/lambda)*(-log(u)) ),
    #                   (1/alpha)*log((alpha/(lambda*exp(beta0))*(-log(u)))- (exp(alpha*v)-1-exp(beta0+alpha*v))/(exp(beta0)))) 
  }
  ##### Exponential Censoring + maximu followup time of 10 years
  cens <- pmin(stats::rexp(n.sample,mu), 10)
  obs.tm <- pmin(event.tm,cens)
  delta <- event.tm < cens
  # obs.tm <- event.tm
  #  delta <- rep(1, n.sample)
  ### Interval censoring of V ### Create questionairre times
  w <- matrix(nrow = n.sample, ncol = n.points)
  points <- seq(0,10, length.out = n.points+1)
  #   intervals <- seq(0,max(event.tm), length.out = n.intervals+1)
  #  intervals <- seq(0,max(event.tm), by = 5)
  for (j in 1:n.points)
  {
    w[,j] <- stats::runif(n = n.sample, min = points[j], max = points[j+1])
    #   w[,j] <- intervals[j+1]
  }
  
  # w.res is the vector of zeroes and ones (X_i(t)) where timings are the values in w
  w.res <- apply(w,2,function(x) x>v)
  w <- apply(w, 2, function(x) ifelse(x < obs.tm, x , Inf)) # Correct questionairre times after event
  w.res[w==Inf] <- Inf # Correct questionairre times after event
  
  list.back <- list(event.tm = event.tm,obs.tm = obs.tm, delta = delta, 
                    w = w, w.res = w.res, v = v)
  
  return(list.back)
}


# @importFrom msm rpexp
SimCoxIntervalCensSinglePexp <- function(n.sample, lambda = 0.25, alpha= NULL, beta0 = log(0.75), 
                                           mu = 1, n.points = 5, rates, ts)
{
  ## Simulate changepoint time
  #v <- rweibull(n.sample, shape = weib.shape, scale = weib.scale)
  v <- msm::rpexp(n.sample, rate=rates, t=ts)
  #v <- ifelse(v < 0, 0, v)
  #v <- ifelse(v > 10, 10, v)
  #v <- runif(n.sample,0,10)
  #####If alpha is null, Create a piecewise exponential where the intervals are [0,V) and [V,\infty] 
  u <- stats::runif(n.sample) # Random uniform variable
  event.tm <- vector(length = n.sample)
  if (is.null(alpha))
  {
    tm.cond <- -log(u) <  lambda*v # the probability of the condition is the probability of event in the second function piece
    event.tm[tm.cond==1] <- -log(u[tm.cond==1])/lambda
    event.tm[tm.cond==0] <- -(log(u[tm.cond==0]) + lambda*v[tm.cond==0]*(1-exp(beta0)))/(lambda*exp(beta0))
    
    #event.tm <- ifelse(tm.cond, -log(u)/lambda,  -(log(u) + lambda*v*(1-exp(beta0)))/(lambda*exp(beta0)))
  } else
    ##### Create a piecewise Gompertz distribution where the intervals are [0,V) and [V,\infty] (See Austin, 2011)
    
  {
    tm.cond <- -log(u)  < (lambda/alpha)*(exp(alpha*v)-1)
    event.tm[tm.cond==1] <- (1/alpha)*log(1 + (alpha/lambda)*(-log(u[tm.cond==1])) )
    event.tm[tm.cond==0] <- (1/alpha)*log((alpha/(lambda*exp(beta0))*(-log(u[tm.cond==0])))- 
                                            (exp(alpha*v[tm.cond==0])-1-exp(beta0+alpha*v[tm.cond==0]))/(exp(beta0))) 
    #event.tm <- ifelse(tm.cond, (1/alpha)*log(1 + (alpha/lambda)*(-log(u)) ),
    #                   (1/alpha)*log((alpha/(lambda*exp(beta0))*(-log(u)))- (exp(alpha*v)-1-exp(beta0+alpha*v))/(exp(beta0)))) 
  }
  ##### Exponential Censoring + maximu followup time of 5 years
  cens <- pmin(stats::rexp(n.sample,mu), 5)
  obs.tm <- pmin(event.tm,cens)
  delta <- event.tm < cens
  # obs.tm <- event.tm
  #  delta <- rep(1, n.sample)
  ### Interval censoring of V ### Create questionairre times
  w <- matrix(nrow = n.sample, ncol = n.points)
  points <- seq(0,5, length.out = n.points+1)
  #   intervals <- seq(0,max(event.tm), length.out = n.intervals+1)
  #  intervals <- seq(0,max(event.tm), by = 5)
  for (j in 1:n.points)
  {
    w[,j] <- stats::runif(n = n.sample, min = points[j], max = points[j+1])
    #   w[,j] <- intervals[j+1]
  }
  
  # w.res is the vector of zeroes and ones (X_i(t)) where timings are the values in w
  w.res <- apply(w,2,function(x) x>v)
  w <- apply(w, 2, function(x) ifelse(x < obs.tm, x , Inf)) # Correct questionairre times after event
  w.res[w==Inf] <- Inf # Correct questionairre times after event
  
  list.back <- list(event.tm = event.tm,obs.tm = obs.tm, delta = delta, 
                    w = w, w.res = w.res, v = v)
  
  return(list.back)
}
### Auxiliry function to calculate (F(t|x)) (Equation 7 in Wng et Al. Biom 2016)
# On June 29 - just fixing notations


AuxF <- function(tt, q1, q2, beta.q1 = log(2), beta.q2 = log(0.5), u)
{
  Lam <- (log(1+tt) + tt^(1/2))/5
  FF <- 1 - exp(-Lam*exp(q1*beta.q1 + q2*beta.q2))
  return(FF - u)
}


SimCoxIntervalCensCox <- function(n.sample, lambda = 0.25, alpha = 0.1  ,beta0 = log(0.75), gamma.q = c(log(0.75), log(2.5)),
                                  gamma.z = log(1.5), mu = 1, n.points = 5)
{
  q1 <- stats::rbinom(n = n.sample, size = 1, 0.5)
  q2 <- stats::rnorm(n = n.sample, mean = 0, sd =  0.5)
  z1 <- stats::rnorm(n = n.sample, mean = 0, sd =  1)
  ## Simulate changepoint time -- Adopting Wang et Al. simulation exapmle ##
  u.for.v <- stats::runif(n.sample)
  v <- vector(length = n.sample) 
  for(i in 1:n.sample)
  {
    while(sign(AuxF(tt = 0, q1 = q1[i], q2 = q2[i], u = u.for.v[i]))==sign(AuxF(tt = 10000, q1 = q1[i], q2 = q2[i], u = u.for.v[i])))
    {
      q1[i] <- q1[i]*0.9
      q2[i] <- q2[i]*0.9
    }
  v[i] <- stats::uniroot(AuxF, q1 = q1[i], q2 = q2[i], u = u.for.v[i], c(0,10000))$root
  }
  #####If alpha is null, Create a piecewise exponential where the intervals are [0,V) and [V,\infty] 
  u <- stats::runif(n.sample) # Random uniform variable
  event.tm <- vector(length = n.sample)
  ##### Create a piecewise Gompertz distribution where the intervals are [0,V) and [V,\infty] (See Austin, 2011)
  exp.qz.gamma <- exp(q1*gamma.q[1] + q2*gamma.q[2] + z1*gamma.z)
  tm.cond <- -log(u)  < (lambda*exp.qz.gamma/alpha)*(exp(alpha*v)-1)
    event.tm[tm.cond==1] <- (1/alpha)*log(1 + (alpha/(exp.qz.gamma[tm.cond==1]*lambda))*(-log(u[tm.cond==1])) )
    event.tm[tm.cond==0] <- (1/alpha)*log((alpha/(lambda*exp(beta0)*exp.qz.gamma[tm.cond==0])*(-log(u[tm.cond==0])))- 
                                            (exp(alpha*v[tm.cond==0])-1-exp(beta0+alpha*v[tm.cond==0]))/(exp(beta0))) 
  ##### Exponential Censoring + maximum followup time of 5 years
  cens <- pmin(stats::rexp(n.sample,mu), 5)
  obs.tm <- pmin(event.tm,cens)
  delta <- event.tm < cens
  ### Interval censoring of V ### Create questionairre times
  w <- matrix(nrow = n.sample, ncol = n.points)
  points <- seq(0,5, length.out = n.points+1)
  for (j in 1:n.points)
  {
    w[,j] <- stats::runif(n = n.sample, min = points[j], max = points[j+1])
    #   w[,j] <- intervals[j+1]
  }
  
  # w.res is the vector of zeroes and ones (X_i(t)) where timings are the values in w
  w.res <- apply(w,2,function(x) x>v)
  w <- apply(w, 2, function(x) ifelse(x < obs.tm, x , Inf)) # Correct questionairre times after event
  w.res[w==Inf] <- Inf # Correct questionairre times after event
  
  list.back <- list(event.tm = event.tm,obs.tm = obs.tm, delta = delta, 
                    w = w, w.res = w.res, v = v, Z = cbind(q1, q2, z1), Q = cbind(q1, q2))
  
  return(list.back)
}
## July 2017 - non-terminal event  to make sure MidI works then under the null ###


SimCoxIntervalCensCoxNonTerminal <- function(n.sample, lambda = 0.25, alpha = 0.1  ,beta0 = log(0.75), gamma.q = c(log(0.75), log(2.5)),
                                  gamma.z = log(1.5), mu = 1, n.points = 5)
{
  q1 <- stats::rbinom(n = n.sample, size = 1, 0.5)
  q2 <- stats::rnorm(n = n.sample, mean = 0, sd =  0.5)
  z1 <- stats::rnorm(n = n.sample, mean = 0, sd =  1)
  ## Simulate changepoint time -- Adopting Wang et Al. simulation exapmle ##
  u.for.v <- stats::runif(n.sample)
  v <- vector(length = n.sample) 
  for(i in 1:n.sample)
  {
    while(sign(AuxF(tt = 0, q1 = q1[i], q2 = q2[i], u = u.for.v[i]))==sign(AuxF(tt = 10000, q1 = q1[i], q2 = q2[i], u = u.for.v[i])))
    {
      q1[i] <- q1[i]*0.9
      q2[i] <- q2[i]*0.9
    }
    v[i] <- stats::uniroot(AuxF, q1 = q1[i], q2 = q2[i], u = u.for.v[i], c(0,10000))$root
  }
  #####If alpha is null, Create a piecewise exponential where the intervals are [0,V) and [V,\infty] 
  u <- stats::runif(n.sample) # Random uniform variable
  event.tm <- vector(length = n.sample)
  ##### Create a piecewise Gompertz distribution where the intervals are [0,V) and [V,\infty] (See Austin, 2011)
  exp.qz.gamma <- exp(q1*gamma.q[1] + q2*gamma.q[2] + z1*gamma.z)
  tm.cond <- -log(u)  < (lambda*exp.qz.gamma/alpha)*(exp(alpha*v)-1)
  event.tm[tm.cond==1] <- (1/alpha)*log(1 + (alpha/(exp.qz.gamma[tm.cond==1]*lambda))*(-log(u[tm.cond==1])) )
  event.tm[tm.cond==0] <- (1/alpha)*log((alpha/(lambda*exp(beta0)*exp.qz.gamma[tm.cond==0])*(-log(u[tm.cond==0])))- 
                                          (exp(alpha*v[tm.cond==0])-1-exp(beta0+alpha*v[tm.cond==0]))/(exp(beta0))) 
  ##### Exponential Censoring + maximum followup time of 5 years
  cens <- pmin(stats::rexp(n.sample,mu), 5)
  obs.tm <- pmin(event.tm,cens)
  delta <- event.tm < cens
  ### Interval censoring of V ### Create questionairre times
  w <- matrix(nrow = n.sample, ncol = n.points)
  points <- seq(0,5, length.out = n.points+1)
  for (j in 1:n.points)
  {
    w[,j] <- stats::runif(n = n.sample, min = points[j], max = points[j+1])
    #   w[,j] <- intervals[j+1]
  }
  
  # w.res is the vector of zeroes and ones (X_i(t)) where timings are the values in w
  w.res <- apply(w,2,function(x) x>v)
  # commented out the next two lines fo non-terminal events
  #w <- apply(w, 2, function(x) ifelse(x < obs.tm, x , Inf)) # Correct questionairre times after event
  #w.res[w==Inf] <- Inf # Correct questionairre times after event
  
  list.back <- list(event.tm = event.tm,obs.tm = obs.tm, delta = delta, 
                    w = w, w.res = w.res, v = v, Z = cbind(q1, q2, z1), Q = cbind(q1, q2))
  
  return(list.back)
}


