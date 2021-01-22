## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(ssdtools)
library(ggplot2)

## ---- message=FALSE-----------------------------------------------------------
library(purrr)
library(tidyr)
library(dplyr)

boron_preds <- nest(ssdtools::boron_data, data = c(Chemical, Species, Conc, Units)) %>%
  mutate(
    Fit = map(data, ssd_fit_dists, dists = "lnorm"),
    Prediction = map(Fit, predict)
  ) %>%
  unnest(Prediction)

## ---- fig.width = 5, fig.height = 5-------------------------------------------
ssd_plot(boron_data, boron_preds, xlab = "Concentration (mg/L)", ci = FALSE) +
  facet_wrap(~Group)

## ---- fig.width = 5, fig.height = 5-------------------------------------------
set.seed(10)
ssd_plot_cf(boron_data)

## ---- fig.width=6, fig.height=6, fig.show='hold'------------------------------
plot(boron_dists)

