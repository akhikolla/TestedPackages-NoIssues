## ---- setseed, echo=FALSE-----------------------------------------------------
set.seed(1)
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)
if("package:GillespieSSA" %in% search()) detach("package:GillespieSSA", unload=TRUE) 

## -----------------------------------------------------------------------------
library(GillespieSSA2)
sim_name <- "SIRS metapopulation model"
patchPopSize <- 500                    # Patch size
U <- 20                                # Number of patches
final_time <- 50                       # Final time

params <- c(
  beta = 0.001,                        # Transmission rate
  gamma = 0.1,                         # Recovery rate
  rho = 0.005,                         # Loss of immunity rate
  epsilon = 0.01,                      # Proportion inter-patch transmissions
  N = patchPopSize                     # Patch population size (constant)
) 

## -----------------------------------------------------------------------------
initial_state <- c(patchPopSize - 1, 1, rep(c(patchPopSize, 0), U - 1))
names(initial_state) <- unlist(lapply(seq_len(U), function(i) paste0(c("S", "I"), i)))

## -----------------------------------------------------------------------------
reactions <- unlist(lapply(
  seq_len(U),
  function(patch) {
    i <- patch
    j <- if (patch == 1) U else patch - 1
    
    Si <- paste0("S", i)
    Ii <- paste0("I", i)
    Ij <- paste0("I", j)
    list(
      reaction(
        propensity = paste0("(1 - epsilon) * beta * ", Si, " * ", Ii), 
        effect = setNames(c(-1, +1), c(Si, Ii)),
        name = paste0("intra_patch_infection_", i)
      ),
      reaction(
        propensity = paste0("epsilon * beta * ", Si, " * ", Ij),
        effect = setNames(c(-1, +1), c(Si, Ii)),
        name = paste0("inter_patch_infection_", i)
      ), 
      reaction(
        propensity = paste0("gamma * ", Ii),
        effect = setNames(-1, Ii),
        name = paste0("recovery_from_infection_", i)
      ),
      reaction(
        propensity = paste0("rho * (N - ", Si, " - ", Ii, ")"),
        effect = setNames(+1, Si),
        name = paste0("loss_of_immunity_", i)
      )
    )
  }
), recursive = FALSE)

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
  method = ssa_etl(),
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

