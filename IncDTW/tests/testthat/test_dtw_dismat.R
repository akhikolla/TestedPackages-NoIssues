context("dtw_dismat")
# library(dtw)

test_that("correct values univariate", {
   # skip("")
   eps <- 1/10^7
   n_min <- 20
   n_max <- 30
   WS <- abs(n_min - n_max) + 1
   set.seed(1210)
   lot <- lapply(1:10, function(i){
      rnorm(sample(seq(n_min, n_max, 2), 1), log(i+1))
   })
   dm1 <- dtw_dismat(lot, ws = WS, dist_method = "norm1", ncores = 1, return_matrix = FALSE)
   dm2 <- dtw_dismat(lot, ws = WS, dist_method = "norm1", ncores = 2, return_matrix = FALSE)
   
   tmp <- sample(10, 2, replace=F)
   i <- tmp[1]
   j <- tmp[2]
   byhand <- sapply(1:(length(lot)-1), function(i){
      sapply((i+1):length(lot), function(j){
         dtw2vec(lot[[i]], lot[[j]], dist_method = "norm1", ws = WS)$normalized_distance
      })
   })
   byhand <- unlist(byhand)
   
   sum(abs(byhand-dm1$dismat))
   sum(abs(byhand-dm2$dismat))
   
   expect_equal(sum(abs(byhand-dm1$dismat)) < eps, TRUE)
   expect_equal(sum(abs(byhand-dm2$dismat)) < eps, TRUE)
   
})



test_that("correct values multivariate, not normalized", {
   # skip("")
   eps <- 1/10^7
   n_min <- 20
   n_max <- 30
   WS <- 1000#NULL#abs(n_min - n_max) + 1
   SP <- "symmetric2"
   DM <- "norm2"
   
   lot <- lapply(1:3, function(i){
      matrix(rnorm(sample(seq(n_min*2, n_max*2, 2), 1), log(i+1)), ncol = 2)
   })
   dm1 <- dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = FALSE,
                     ncores = 1, return_matrix = FALSE)
   dm2 <- dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = FALSE,
                     ncores = 2, return_matrix = FALSE)
   
   # tmp <- sample(10, 2, replace=F)
   # i <- tmp[1]
   # j <- tmp[2]
   byhand <- sapply(1:(length(lot)-1), function(i){
      sapply((i+1):length(lot), function(j){
         dtw2vec(lot[[i]], lot[[j]], step_pattern = SP, dist_method = DM, ws = WS)$distance
      })
   })
   byhand <- unlist(byhand)
   
   sum(abs(byhand-dm1$dismat))
   sum(abs(byhand-dm2$dismat))
   
   expect_equal(sum(abs(byhand-dm1$dismat)) < eps, TRUE)
   expect_equal(sum(abs(byhand-dm2$dismat)) < eps, TRUE)
   
})



test_that("correct values multivariate, normalized", {
   # skip("")
   eps <- 1/10^7
   n_min <- 20
   n_max <- 30
   WS <- 1000#NULL#abs(n_min - n_max) + 1
   SP <- "symmetric2"
   DM <- "norm2"
   
   lot <- lapply(1:3, function(i){
      matrix(rnorm(sample(seq(n_min*2, n_max*2, 2), 1), log(i+1)), ncol = 2)
   })
   dm1 <- dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = TRUE,
                     ncores = 1, return_matrix = FALSE)
   dm2 <- dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = TRUE,
                     ncores = 2, return_matrix = FALSE)
   
   # tmp <- sample(10, 2, replace=F)
   # i <- tmp[1]
   # j <- tmp[2]
   byhand <- sapply(1:(length(lot)-1), function(i){
      sapply((i+1):length(lot), function(j){
         dtw2vec(lot[[i]], lot[[j]], step_pattern = SP, dist_method = DM, ws = WS)$normali
      })
   })
   byhand <- unlist(byhand)
   
   sum(abs(byhand-dm1$dismat))
   sum(abs(byhand-dm2$dismat))
   
   expect_equal(sum(abs(byhand-dm1$dismat)) < eps, TRUE)
   expect_equal(sum(abs(byhand-dm2$dismat)) < eps, TRUE)
   
})



test_that("for dev only, speed comparison", {
   skip("")
   eps <- 1/10^7
   n_min <- 20
   n_max <- 30
   WS <- 1000#NULL#abs(n_min - n_max) + 1
   SP <- "symmetric1"
   DM <- "norm2"
   NORM <- FALSE
   
   lot <- lapply(1:3, function(i){
      # matrix(rnorm(40, log(i+1)), ncol = 2)
      matrix(rnorm(sample(seq(n_min*2, n_max*2, 2), 1), log(i+1)), ncol = 2)
   })
   lot_t <- lapply(lot, t)
   
   f0 <- function(lot){dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = NORM,
                     ncores = 1, return_matrix = FALSE)$dismat}
   fcp <- function(lot){dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = NORM,
                     ncores = 3, return_matrix = FALSE, useRcppParallel = TRUE)$dismat}
   frp <- function(lot){dtw_dismat(lot, ws = WS, step_pattern = SP, dist_method = DM, normalize = NORM,
                     ncores = 3, return_matrix = FALSE, useRcppParallel = FALSE)$dismat}
   fpp <- function(lot){parallelDist::parDist(x=lot, ws = WS, method = "dtw", 
                                              step_pattern = dtw::symmetric1, normalize = "blabla",threads = 3)}
   
   f0(lot)
   fcp(lot)
   frp(lot)
   fpp(lot_t)
   # conclusion: parallelDist does always use step_pattern == symmetric1, try it out
   
   tic <- Sys.time()
   nlot <- 100# number of time series
   tmp <- lapply(seq(250, 1000, 250), function(nn){
      set.seed(nn)
      
      lot <- lapply(1:nlot, function(i){
         matrix(rnorm(nn*2, log(i),1), ncol = 2)
      })
      lot_t <- lapply(lot, t)
      
      mic <- microbenchmark::microbenchmark(f0(lot),
                                            fcp(lot),
                                            # frp(lot),
                                            fpp(lot_t),
                                            times = 2
      )
      
      return(data.frame(expr = mic$expr, times = mic$time, nn = nn))
   })
   tac <- Sys.time()
   tac-tic
   mics <- do.call(rbind, tmp)
   class(mics) <- c("microbenchmark", "data.frame")
   dfp <- aggregate(times ~ nn+expr, data = mics, median)
   dfp$sec <- dfp$times/10^9
   
   require(ggplot2)
   ggplot(dfp)+geom_line(aes(x=nn, y = sec, group = expr, col=expr))+scale_y_log10()
   require(data.table)
   dfp2 <- as.data.table(dfp)
   dfp2 <- dfp2[, list(expr, ratio = times/.SD[expr == "fcp(lot)", times]), by = c("nn") ]
   ggplot(dfp2)+geom_line(aes(x=nn, y = ratio, group = expr, col=expr))
})




