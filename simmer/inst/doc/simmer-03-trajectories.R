## ---- cache = FALSE, include=FALSE--------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>", 
                      fig.width = 6, fig.height = 4, fig.align = "center")

## ---- message=FALSE-----------------------------------------------------------
library(simmer)

## -----------------------------------------------------------------------------
traj <- trajectory() %>%
  seize(resource = "doctor", amount = 1) %>%
  timeout(task = 3) %>%
  release(resource = "doctor", amount = 1)

## -----------------------------------------------------------------------------
methods(class="trajectory")

## -----------------------------------------------------------------------------
trajectory() %>%
  timeout(rexp(1, 10)) %>%        # fixed
  timeout(function() rexp(1, 10)) # dynamic

## ---- error = TRUE------------------------------------------------------------
traj <- trajectory() %>%
  log_(function() as.character(now(env)))

env <- simmer() %>%
  add_generator("dummy", traj, function() 1) %>%
  run(4)

## -----------------------------------------------------------------------------
traj <- trajectory() %>%
  log_(function() as.character(now(env)))

env <- simmer() %>%
  add_generator("dummy", traj, function() 1)

env %>% run(4) %>% invisible

## -----------------------------------------------------------------------------
# first, instantiate the environment
env <- simmer()

# here I'm using it
traj <- trajectory() %>%
  log_(function() as.character(now(env)))

# and finally, run it
env %>%
  add_generator("dummy", traj, function() 1) %>%
  run(4) %>% invisible

## -----------------------------------------------------------------------------
t1 <- trajectory() %>% seize("dummy", 1)
t2 <- trajectory() %>% timeout(1)
t3 <- trajectory() %>% release("dummy", 1)

t0 <- join(t1, t2, t3)
t0

## -----------------------------------------------------------------------------
t0 <- trajectory() %>%
  join(t1) %>%
  timeout(1) %>%
  join(t3)
t0

## -----------------------------------------------------------------------------
length(t0)

## -----------------------------------------------------------------------------
t0[c(TRUE, FALSE, TRUE)]

## -----------------------------------------------------------------------------
t0[c(1, 3)]
t0[c(3, 1)]

## -----------------------------------------------------------------------------
t0[-2]

## -----------------------------------------------------------------------------
t0[c("seize", "release")]
t0[c("release", "seize")]

## -----------------------------------------------------------------------------
t0[]

## -----------------------------------------------------------------------------
head(t0, 2)
tail(t0, -1)

## -----------------------------------------------------------------------------
t0[[2]]

## -----------------------------------------------------------------------------
join(t0, t0)["timeout"]
join(t0, t0)[["timeout"]]

