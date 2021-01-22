## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo=T,cache=T,results='hold');
library(Temporal);

## -----------------------------------------------------------------------------
# Generate exponential time to event data
data = genData(n=1e3,dist="exp",theta=c(2),p=0.2)
# Estimate parameters
fit= fitParaSurv(time=data$time,status=data$status,dist="exp");
show(fit);

## -----------------------------------------------------------------------------
# Generate gamma time to event data
data = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.25);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="gamma",tau=0.5);
show(fit);

## -----------------------------------------------------------------------------
set.seed(102);
# Generate generalized gamma time to event data
data = genData(n=1e4,dist="gen-gamma",theta=c(2,2,2),p=0.1);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="gen-gamma",report=T);
show(fit);

## ---- results='markup', eval=FALSE--------------------------------------------
#  # Initialization
#  fitParaSurv(time=data$time,status=data$status,dist="gen-gamma",init=c(2,2,2));

## -----------------------------------------------------------------------------
# Generate log-normal time to event data
data = genData(n=1e3,dist="log-normal",theta=c(1,2),p=0.15);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="log-normal",tau=c(5,10,25));
show(fit);

## -----------------------------------------------------------------------------
# Generate Weibull time to event data
data = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.3);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="weibull");
show(fit);

## -----------------------------------------------------------------------------
set.seed(101);
# Target group
d1 = genData(n=1e3,dist="gamma",theta=c(2,1),p=0.25);
d1$arm = 1;
# Reference group 
d0 = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.15);
d0$arm = 0;
# Overall data 
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="gamma",dist0="gamma");
cat("\n");
show(comp);

## -----------------------------------------------------------------------------
# Target group
d1 =genData(n=1e3,dist="weibull",theta=c(2,2),p=0.5);
d1$arm = 1;
# Reference group
d0 = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.0);
d0$arm = 0;
# Overall data
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="weibull",dist0="weibull",boot=1e3);
cat("\n");
show(comp);

## -----------------------------------------------------------------------------
set.seed(105);
# Target group
d1 = genData(n=1e3,dist="log-normal",theta=c(0,sqrt(2*log(2))),p=0.1);
d1$arm = 1;
# Reference group
d0 = genData(n=1e3,dist="weibull",theta=c(2,sqrt(log(2))),p=0.1);
d0$arm = 0;
# Overall data
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="log-normal",dist0="weibull",perm=1e3);
cat("\n");
show(comp);

## -----------------------------------------------------------------------------
set.seed(106);
# Target group
d1 = genData(n=1e3,dist="gamma",theta=c(4,4),p=0.2);
d1$arm = 1;
# Reference group
d0 = genData(n=1e3,dist="exp",theta=c(1),p=0.2);
d0$arm = 0;
# Overall data
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="gamma",dist0="exp",tau=c(0.5,1.0,1.5),perm=1e3);
cat("\n");
show(comp);

