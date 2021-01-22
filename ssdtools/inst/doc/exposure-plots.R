## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.width = 6,
  fig.height = 4
)

## ----echo=TRUE, warning=FALSE, message=FALSE----------------------------------
library(ggplot2)
library(ssdtools)
data(boron_data)
fit <- ssd_fit_dists(boron_data, dists = c("llogis", "lnorm", "gamma"))
fit.plot <- autoplot(fit)
fit.plot

## ----echo=FALSE---------------------------------------------------------------
# what are the mean and sd on the log-scale for the envirionmental concentration
ex.mean.log <- log(0.1)
ex.sd.log <- 1

## -----------------------------------------------------------------------------
ex.cdf <- data.frame(Conc = exp(seq(log(.01), log(10), .1))) # generate a grid of concentrations
ex.cdf$ex.cdf <- plnorm(ex.cdf$Conc,
  meanlog = ex.mean.log,
  sdlog = ex.sd.log
) # generate the cdf

## -----------------------------------------------------------------------------
fit.plot + geom_line(data = ex.cdf, aes(x = Conc, y = ex.cdf), color = "red", size = 2) +
  annotate("text",
    label = paste("Exposure distribution"),
    x = 1.08 * ex.cdf$Conc[which.max(ex.cdf$ex.cdf > 0.5)], y = 0.5, angle = 75
  )

## ----echo=TRUE----------------------------------------------------------------
set.seed(99)
ex.risk <- ssd_exposure(fit, meanlog = ex.mean.log, sdlog = ex.sd.log)
ex.risk

## ----echo=TRUE----------------------------------------------------------------
fit.plot + geom_line(dat = ex.cdf, aes(x = Conc, y = ex.cdf), color = "red", size = 2) +
  annotate("text",
    label = paste("Exposure distribution"),
    x = 1.08 * ex.cdf$Conc[which.max(ex.cdf$ex.cdf > 0.5)], y = 0.5, angle = 75
  ) +
  annotate("text",
    label = paste("Verdonck risk :", round(ex.risk, 5)),
    x = Inf, y = 0, hjust = 1.1, vjust = -.5
  )

