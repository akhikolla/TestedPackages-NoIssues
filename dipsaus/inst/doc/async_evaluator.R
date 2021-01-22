## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = FALSE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#  # Future
#  plan(multisession, workers = 2L)
#  
#  start = Sys.time()
#  lapply(1:8, function(ii){
#    future({ Sys.sleep(2) })
#    print(sprintf('%d - Time ellapsed: %.2f sec.',
#                  ii, time_delta(start, Sys.time())))
#  })
#  #> [1] "1 - Time ellapsed: 0.22 sec."
#  #> [1] "2 - Time ellapsed: 0.43 sec."
#  #> [1] "3 - Time ellapsed: 2.47 sec."
#  #> [1] "4 - Time ellapsed: 2.67 sec."
#  #> [1] "5 - Time ellapsed: 4.72 sec."
#  #> [1] "6 - Time ellapsed: 4.94 sec."
#  #> [1] "7 - Time ellapsed: 6.98 sec."
#  #> [1] "8 - Time ellapsed: 7.19 sec."

## -----------------------------------------------------------------------------
#  evaluator <- make_async_evaluator(name = 'test', n_nodes = 1, n_subnodes = 7)
#  #> ✔ Initializing 2 workers with 1 sub-workers. Total workers will be (0+2) x 1

## -----------------------------------------------------------------------------
#  evaluator$run({ Sys.getpid() })

## -----------------------------------------------------------------------------
#  evaluator$progress()
#  #> [1]  1  0  0  1

## -----------------------------------------------------------------------------
#  evaluator$run({ Sys.getpid() }, success = print)
#  #> [1] 21304

## -----------------------------------------------------------------------------
#  start = Sys.time()
#  system.time({
#    lapply(1:20, function(ii){
#      evaluator$run({
#        Sys.sleep(10)
#      },success = function(v){
#        print(Sys.time() - start)
#      })
#    })
#  })
#  #>    user  system elapsed
#  #>   0.073   0.020   0.097
#  
#  ## You can still run anything within main session
#  ## Print current progress
#  ## Total  Running  await  finished
#  print(evaluator$progress())
#  #> [1] 22  1  21  0
#  #> Time difference of 13.11858 secs
#  #> Time difference of 13.11977 secs
#  #> Time difference of 13.12027 secs
#  #> Time difference of 13.12076 secs
#  #> Time difference of 13.12104 secs
#  #> Time difference of 13.12128 secs
#  #> Time difference of 14.30096 secs
#  #> Time difference of 22.93647 secs
#  #> Time difference of 22.93681 secs
#  #> Time difference of 24.45621 secs
#  #> Time difference of 24.45656 secs
#  #> Time difference of 24.45683 secs
#  #> Time difference of 24.45708 secs
#  #> Time difference of 25.7501 secs
#  #> Time difference of 33.89107 secs
#  #> Time difference of 34.90931 secs
#  #> Time difference of 34.90963 secs
#  #> Time difference of 35.93507 secs
#  #> Time difference of 35.9354 secs
#  #> Time difference of 35.93566 secs

## -----------------------------------------------------------------------------
#  pid1 <- Sys.getpid()
#  evaluator$run({
#    sprintf('Finished in %s, scheduled in master session %s',
#            Sys.getpid(), !!pid1)
#  },success = print)

## -----------------------------------------------------------------------------
#  ## 800 MB big array
#  dat <- array(rnorm(100000000), c(100,100,100,10,10))

## -----------------------------------------------------------------------------
#  ## [Bad]
#  evaluator$run({
#    apply(!!dat, 5, median)
#  },success = print)
#  #> Error in checkForRemoteErrors(lapply(cl, recvResult)) :
#  #> ...

## -----------------------------------------------------------------------------
#  ## [OK] Pass parameters to ..., or .list
#  start = Sys.time()
#  evaluator$run({
#    apply(dat, 5, median)
#  },success = function(v){
#    print(sprintf('Total run time: %s', Sys.time() - start))
#  }, dat = dat)
#  print(sprintf('Schedule time: %s', Sys.time() - start))
#  #> [1] "Schedule time: 10.5080621242523"
#  #> [1] "Total run time: 18.8798191547394"

## -----------------------------------------------------------------------------
#  e <- make_async_evaluator(name = 'test')
#  identical(e$get_instance(), evaluator$get_instance())

## -----------------------------------------------------------------------------
#  ## Scale down to 1 manager, 2 workers
#  evaluator$scale_down(n_nodes = 1, n_subnodes = 2)
#  #> [1] 1
#  
#  ## Scale up to 1 manager, 7 workers
#  evaluator$scale_up(n_nodes = 1, n_subnodes = 7)
#  #> [1] 1

## -----------------------------------------------------------------------------
#  evaluator$suspend()

## -----------------------------------------------------------------------------
#  evaluator$stop()

## -----------------------------------------------------------------------------
#  evaluator$terminate()

