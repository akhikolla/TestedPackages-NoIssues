## ---- setseed, echo=FALSE-----------------------------------------------------
set.seed(1)
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)
if("package:GillespieSSA2" %in% search()) detach("package:GillespieSSA2", unload=TRUE) 

## ----gssa1--------------------------------------------------------------------
library(GillespieSSA)
parms <- c(c1 = 10, c2 = .01, c3 = 10)
tf <- 2                                        # Final time
simName <- "Lotka predator-prey model"         # Name
x0 <- c(Y1=1000, Y2=1000)
nu <- matrix(c(+1, -1, 0, 0, 1, -1), nrow = 2, byrow = TRUE)
a  <- c("c1*Y1", "c2*Y1*Y2","c3*Y2")
out <- GillespieSSA::ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = ssa.d(),
  simName = simName,
  censusInterval = .001,
  verbose = FALSE,
  consoleInterval = 1
) 
ssa.plot(out, show.title = TRUE, show.legend = FALSE)

## ---- setseed2, echo=FALSE----------------------------------------------------
set.seed(1)
if("package:GillespieSSA" %in% search()) detach("package:GillespieSSA", unload=TRUE) 

## ----gssa2--------------------------------------------------------------------
library(GillespieSSA2)
sim_name <- "Lotka Predator-Prey model"
params <- c(c1 = 10, c2 = .01, c3 = 10)
final_time <- 2
initial_state <- c(Y1 = 1000, Y2 = 1000)
reactions <- list(
  reaction("c1 * Y1",      c(Y1 = +1)),
  reaction("c2 * Y1 * Y2", c(Y1 = -1, Y2 = +1)),
  reaction("c3 * Y2",      c(Y2 = -1))
)
out <- GillespieSSA2::ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_exact(),
  census_interval = .001,
  verbose = FALSE,
  sim_name = sim_name
) 
plot_ssa(out)

## ----port---------------------------------------------------------------------
out <- 
  GillespieSSA2::ssa(
    initial_state = x0,
    reactions = port_reactions(x0 = x0, a = a, nu = nu),
    params = parms,
    method = ssa_exact(),
    final_time = tf,
    census_interval = .001,
    verbose = FALSE,
    sim_name = simName
  )
print(out$stats)
plot_ssa(out)

