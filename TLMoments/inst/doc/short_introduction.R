## -----------------------------------------------------------------------------
library(TLMoments)
sessionInfo()

## -----------------------------------------------------------------------------
xvec <- rgev(100, loc = 10, scale = 5, shape = .2)

## -----------------------------------------------------------------------------
TLMoments(xvec)
TLMoments(xvec, leftrim = 0, rightrim = 1, max.order = 2)

## -----------------------------------------------------------------------------
tlm <- TLMoments(xvec)
parameters(tlm, distr = "gev")

tlm <- TLMoments(xvec, rightrim = 1)
parameters(tlm, distr = "gev")

## -----------------------------------------------------------------------------
tlm <- TLMoments(xvec)
quantiles(parameters(tlm, distr = "gev"), c(.9, .99, .999))

tlm <- TLMoments(xvec, rightrim = 1)
quantiles(parameters(tlm, distr = "gev"), c(.9, .99, .999))

## -----------------------------------------------------------------------------
tlm <- TLMoments(rgev(100), leftrim = 0, rightrim = 1)

summary(tlm)
summary(parameters(tlm, "gev"))
summary(quantiles(parameters(tlm, "gev"), .99))

## -----------------------------------------------------------------------------
summary(tlm, ci.level = .95, select = 3:4)
summary(parameters(tlm, "gev"), select = "shape")

## -----------------------------------------------------------------------------
library(magrittr)

TLMoments(xvec, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999)) %>% 
  summary()

## -----------------------------------------------------------------------------
xmat <- matrix(rgev(100), nc = 4)
xvec <- xmat[, 3]
xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
xdat <- data.frame(station = rep(1:4, each = 25), hq = as.vector(xmat))

## -----------------------------------------------------------------------------
TLMoments(xvec, leftrim = 0, rightrim = 1)
TLMoments(xmat, leftrim = 0, rightrim = 1)
TLMoments(xlist, leftrim = 0, rightrim = 1)
TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1)

## -----------------------------------------------------------------------------
TLMoments(xvec, leftrim = 0, rightrim = 1) %>% 
  parameters("gev")
TLMoments(xmat, leftrim = 0, rightrim = 1) %>% 
  parameters("gev")
TLMoments(xlist, leftrim = 0, rightrim = 1) %>% 
  parameters("gev")
TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1) %>% 
  parameters("gev")

TLMoments(xvec, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))
TLMoments(xmat, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>%
  quantiles(c(.99, .999))
TLMoments(xlist, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))
TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))

## -----------------------------------------------------------------------------
data <- matrix(rgev(25 * 10, shape = .2), ncol = 10)
plot(TLMoments(data))
plot(TLMoments(data)) + ggplot2::theme_minimal()
plot(TLMoments(data, rightrim = 1), distr = c("gev", "gpd", "exp", "gum"))

## -----------------------------------------------------------------------------
(tlm <- as.TLMoments(c(14.1, 4.3, 1.32)))
parameters(tlm, distr = "gev")
quantiles(parameters(tlm, distr = "gev"), c(.9, .99, .999))

(param <- as.parameters(loc = 10, scale = 5, shape = .2, distr = "gev"))
quantiles(param, c(.9, .99, .999))
TLMoments(param)
TLMoments(param, rightrim = 1)

