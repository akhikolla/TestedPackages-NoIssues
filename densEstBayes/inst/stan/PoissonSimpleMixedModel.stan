data
      {
         int<lower=1> n;        int<lower=1> ncZ;
         int<lower=0> y[n];     matrix[n,2] X;
         matrix[n,ncZ] Z;       real<lower=0> sigmabeta;
         real<lower=0> ssigma;        
      }
      parameters 
      {
         vector[2] beta;          vector[ncZ] u;
         real<lower=0> sigma;
      }
      transformed parameters
      {
         vector[n] mu;
         mu = exp(X*beta + Z*u);
      }
      model 
      {
         for (i in 1:n)
            y[i] ~ poisson(mu[i]);
         beta ~ normal(0,sigmabeta); 
         u ~ normal(0,sigma); 
         sigma ~ cauchy(0,ssigma);
      }
