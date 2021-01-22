## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 4,
  fig.width = 6
)

## ---- message=FALSE-----------------------------------------------------------
library(ssdtools)
library(ggplot2)

set.seed(7)
ssd_plot_cdf(ssd_match_moments(meanlog = 2, sdlog = 2))

## ----fig.height=5-------------------------------------------------------------
set.seed(7)
ssd_plot_cdf(ssd_match_moments(dists = c(
  "gamma",
  "gompertz", "lgumbel", "llogis",
  "lnorm", "weibull"
), meanlog = 2, sdlog = 2)) +
  theme(legend.position = "bottom")

