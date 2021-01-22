## ---- setseed, echo=FALSE-----------------------------------------------------
set.seed(1)
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)
if("package:GillespieSSA" %in% search()) detach("package:GillespieSSA", unload=TRUE) 

## -----------------------------------------------------------------------------
library(GillespieSSA2)
sim_name <- "Rosenzweig-MacArthur Predator-Prey model"
params <- c(
  b = 2, 
  d = 1,
  K = 1000,
  alpha = 0.005, 
  w = 0.0025,
  c = 2,
  g = 2
)
final_time <- 10
initial_state  <- c(N = 500, P = 500)

## -----------------------------------------------------------------------------
reactions <- list(
  reaction("b * N", c(N = +1)),
  reaction("(d + (b - d) * N / K) * N", c(N = -1)),
  reaction("alpha / (1 + w * N) * N * P", c(N = -1)),
  reaction("c * alpha / ( 1 + w * N) * N * P", c(P = +1)),
  reaction("g * P", c(P = -1))
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
  method = ssa_etl(tau = .01),
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
  method = ssa_btl(),
  sim_name = sim_name
) 
plot_ssa(out)

