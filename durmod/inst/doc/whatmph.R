## ----data---------------------------------------------------------------------
library(durmod)
data(durdata)
head(durdata,15)

## -----------------------------------------------------------------------------
risksets <- list(unemp=c('job','program'), onprogram='job')

## -----------------------------------------------------------------------------
ctrl <- mphcrm.control(iters=4, threads=1)

## -----------------------------------------------------------------------------
set.seed(42) # for reproducibility
opt <- mphcrm(d ~ x1 + x2 + 
                C(job, alpha) + ID(id) + D(duration) + S(state), 
              data=durdata, risksets=risksets, control=ctrl)

## -----------------------------------------------------------------------------
print(opt)

## -----------------------------------------------------------------------------
best <- opt[[1]]
summary(best)

## -----------------------------------------------------------------------------
t(sapply(opt, function(o) summary(o)$coefs["job.alpha",]))

## -----------------------------------------------------------------------------
summary(fit[[1]])

## ----eval=FALSE---------------------------------------------------------------
#  library(durmod)
#  data(durdata)
#  newfit <- eval(attr(fit,'call'))

## -----------------------------------------------------------------------------
round(mphdist(fit[[1]]),6)
# and the moments,
mphmoments(fit[[1]])
# and covariance matrix
mphcov(fit[[1]])
# and some pseudo R2
pseudoR2(fit)

## -----------------------------------------------------------------------------
attributes(durdata)[c('means','cov')]

## -----------------------------------------------------------------------------
summary(fit[[which.min(sapply(fit,AIC))]])

## ----echo=FALSE---------------------------------------------------------------
def <- mphcrm.control()

