## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
#options(tibble.print_min = 4L, tibble.print_max = 4L)
library(ztpln)
set.seed(123)

## ---- eval = T----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(123)

rztpln(n = 10, mu = 0, sig = 1)

rztpln(n = 10, mu = 6, sig = 4)

## ---- eval = T----------------------------------------------------------------
rztplnm(n = 100, 
  mu = c(0, 4), 
  sig = c(0.5, 0.5), 
  theta = c(0.2, 0.8)) %>%
  log %>% hist(main = "", xlab = "k")

## ---- eval = T----------------------------------------------------------------
y <- rztpln(n = 100, mu = 2, sig = 5, type1 = TRUE)
sum(dztpln(y, mu = 2, sig = 5, log = TRUE, type1 = TRUE))
sum(dztpln(y, mu = 2, sig = 5, log = TRUE, type1 = FALSE))

## ---- eval = T----------------------------------------------------------------
y2 <- rztpln(n = 100, mu = 2, sig = 5, type1 = FALSE)
sum(dztpln(y2, mu = 2, sig = 5, log = TRUE, type1 = TRUE))
sum(dztpln(y2, mu = 2, sig = 5, log = TRUE, type1 = FALSE))


## -----------------------------------------------------------------------------
k1 <- 1
k <- 1000
dat <- tibble(type1 = dztpln(k1:k, mu = 1, sig = 2, type1 = TRUE),
         type2 = dztpln(k1:k, mu = 1, sig = 2, type1 = FALSE),
         x = k1:k) %>%
  pivot_longer(-x, names_to = "type", values_to = "y")

ggplot(dat %>% dplyr::filter(x <= 10), aes(x = x, y = y, col = type)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  xlab("k") +
  ylab("Likelihood") +
  theme_bw()

## -----------------------------------------------------------------------------
ggplot(dat, aes(x = x, y = y, col = type)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  xlab("k") +
  ylab("Likelihood") +
  theme_bw()

