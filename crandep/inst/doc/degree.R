## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(crandep)
library(igraph)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
g0.imports <- get_graph_all_packages(type = "imports")
d0.imports <- g0.imports %>% igraph::degree(mode = "in")
df0.imports <-
    data.frame(name = names(d0.imports), degree = as.integer(d0.imports)) %>%
    dplyr::arrange(dplyr::desc(degree), name)
head(df0.imports, 10)

## -----------------------------------------------------------------------------
g0.rev_imports <- get_graph_all_packages(type = "reverse imports")
d0.rev_imports <- g0.rev_imports %>% igraph::degree(mode = "out") # note the difference to above
df0.rev_imports <-
    data.frame(name = names(d0.rev_imports), degree = as.integer(d0.rev_imports)) %>%
    dplyr::arrange(dplyr::desc(degree), name)
head(df0.rev_imports, 10)

## -----------------------------------------------------------------------------
identical(df0.imports, df0.rev_imports)
setdiff(df0.imports, df0.rev_imports)
setdiff(df0.rev_imports, df0.imports)

## -----------------------------------------------------------------------------
df1.imports <- df0.imports %>%
	dplyr::filter(degree > 0L) %>% # to prevent warning when plotting on log-log scale
    dplyr::count(degree, name = "frequency") %>%
	dplyr::arrange(dplyr::desc(degree)) %>%
	dplyr::mutate(survival = cumsum(frequency)/sum(frequency))

## -----------------------------------------------------------------------------
gg0 <- df1.imports %>%
	ggplot2::ggplot() +
	ggplot2::geom_point(aes(degree, frequency), size = 0.75) +
	ggplot2::scale_x_log10() +
	ggplot2::scale_y_log10() +
	ggplot2::coord_cartesian(ylim = c(1L, 1e+3L)) +
	ggplot2::theme_bw(12)
gg0

## ----results = FALSE----------------------------------------------------------
x <- dplyr::filter(df0.imports, degree > 0L)$degree # data
u <- 1L # threshold
xi1 <- 1.0 # initial value
a_xi1 <- 0.0 # lower bound of uniform distribution
b_xi1 <- 100.0 # upper bound of uniform distribution
set.seed(3075L)
mcmc0.imports <- mcmc_upp(x = x, u = u, xi1 = xi1, a_xi1 = a_xi1, b_xi1 = b_xi1) # takes seconds

## -----------------------------------------------------------------------------
mcmc0.imports %>%
    ggplot2::ggplot() +
	ggplot2::geom_density(aes(xi1)) +
	ggplot2::theme_bw(12)

## -----------------------------------------------------------------------------
mcmc0.imports <- mcmc0.imports %>%
    dplyr::mutate(alpha = 1.0 / xi1 + 1.0)
mcmc0.imports %>%
	ggplot2::ggplot() +
	ggplot2::geom_density(aes(alpha)) +
	ggplot2::theme_bw(12)

## -----------------------------------------------------------------------------
n0 <- sum(df1.imports$frequency) # TOTAL number of data points in x
## or n0 <- length(x)
n1 <- length(df1.imports$frequency) # number of UNIQUE data points
N <- length(mcmc0.imports$xi1)
freq0 <- surv0 <- matrix(as.numeric(NA), N, n1)
for (i in seq(N)) {
    freq0[i,] <- n0 * dupp(x = df1.imports$degree, u = 1L, xi1 = mcmc0.imports$xi1[i])
    surv0[i,] <- Supp(x = df1.imports$degree, u = 1L, xi1 = mcmc0.imports$xi1[i])
}
df1.imports <- df1.imports %>%
    dplyr::mutate(
	    frequency.mean = apply(freq0, 2, mean),
		frequency.qlow = apply(freq0, 2, quantile, p = 0.025),
		frequency.qupp = apply(freq0, 2, quantile, p = 0.975),
		survival.mean = apply(surv0, 2, mean),
		survival.qlow = apply(surv0, 2, quantile, p = 0.025),
		survival.qupp = apply(surv0, 2, quantile, p = 0.975)
	)

## -----------------------------------------------------------------------------
gg1 <- df1.imports %>%
    ggplot2::ggplot() +
	ggplot2::geom_point(aes(degree, frequency), size = 0.75) +
	ggplot2::geom_line(aes(degree, frequency.mean), col = 4, lty = 2) +
	ggplot2::geom_line(aes(degree, frequency.qlow), col = 2, lty = 3) +
	ggplot2::geom_line(aes(degree, frequency.qupp), col = 2, lty = 3) +
	ggplot2::scale_x_log10() +
	ggplot2::scale_y_log10() +
	ggplot2::coord_cartesian(ylim = c(1L, 1e+3L)) +
	ggplot2::theme_bw(12)
gg1

## -----------------------------------------------------------------------------
gg2 <- df1.imports %>%
    ggplot2::ggplot() +
	ggplot2::geom_point(aes(degree, survival), size = 0.75) +
	ggplot2::geom_line(aes(degree, survival.mean), col = 4, lty = 2) +
	ggplot2::geom_line(aes(degree, survival.qlow), col = 2, lty = 3) +
	ggplot2::geom_line(aes(degree, survival.qupp), col = 2, lty = 3) +
	ggplot2::scale_x_log10() +
	ggplot2::scale_y_log10() +
	ggplot2::theme_bw(12)
gg2

