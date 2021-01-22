## ---- setseed, echo=FALSE-----------------------------------------------------
set.seed(1)
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)
if("package:GillespieSSA" %in% search()) detach("package:GillespieSSA", unload=TRUE) 

## -----------------------------------------------------------------------------
library(GillespieSSA2)
initial_state <- c(prey = 1000, predators = 1000)

## -----------------------------------------------------------------------------
params <- c(c1 = 10, c2 = 0.01, c3 = 10)
reactions <- list(
  #        propensity function        effects                          name for reaction
  reaction(~c1 * prey,                c(prey = +1),                    name = "prey_up"),
  reaction(~c2 * prey * predators,    c(prey = -1, predators = +1),    name = "predation"),
  reaction(~c3 * predators,           c(predators = -1),               name = "pred_down")
)

## -----------------------------------------------------------------------------
out <- 
  ssa(
    initial_state = initial_state,
    reactions = reactions,
    params = params,
    method = ssa_exact(),
    final_time = 5,
    census_interval = .001,
    verbose = TRUE,
    sim_name = "My first SSA run!"
  )

## -----------------------------------------------------------------------------
print(out$stats)

## ----fig.width=7, fig.height=5------------------------------------------------
plot_ssa(out)

