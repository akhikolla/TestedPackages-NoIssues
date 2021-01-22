## ---- setseed, echo=FALSE-----------------------------------------------------
set.seed(1)
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)
if("package:GillespieSSA" %in% search()) detach("package:GillespieSSA", unload=TRUE) 

## -----------------------------------------------------------------------------
library(GillespieSSA2)
sim_name <- "Linear Chain System"
M <- 50
params <- c(c = 1)
final_time <- 5
initial_state <- c(1000, rep(0, M)) 
names(initial_state) <- paste0("x", seq_len(M+1))

## -----------------------------------------------------------------------------
reactions <- lapply(
  seq_len(M),
  function(i) {
    effect <- c(-1, 1)
    names(effect) <- paste0("x", c(i, i + 1))
    
    reaction(paste0("c * x", i), effect)
  }
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
  method = ssa_etl(tau = .1),
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
  method = ssa_btl(mean_firings = 50),
  sim_name = sim_name
) 
plot_ssa(out)

