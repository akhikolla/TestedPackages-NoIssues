## ---- setseed, echo=FALSE-----------------------------------------------------
set.seed(1)
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)
if("package:GillespieSSA" %in% search()) detach("package:GillespieSSA", unload=TRUE) 

## -----------------------------------------------------------------------------
library(GillespieSSA2)
sim_name <- "Lotka Predator-Prey model"
params <- c(c1 = 10, c2 = .01, c3 = 10)
final_time <- 2
initial_state <- c(Y1 = 1000, Y2 = 1000)

## -----------------------------------------------------------------------------
reactions <- list(
  reaction("c1 * Y1", c(Y1 = +1)),
  reaction("c2 * Y1 * Y2", c(Y1 = -1, Y2 = +1)),
  reaction("c3 * Y2", c(Y2 = -1))
)

## ----exact--------------------------------------------------------------------
set.seed(1)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_exact(),
  sim_name = sim_name
) 
plot_ssa(out)

## ----etl----------------------------------------------------------------------
set.seed(1)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_etl(tau = .002),
  sim_name = sim_name
) 
plot_ssa(out)

## ----btl----------------------------------------------------------------------
set.seed(1)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_btl(mean_firings = 100),
  sim_name = sim_name
) 
plot_ssa(out)

