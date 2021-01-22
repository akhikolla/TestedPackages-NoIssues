context("Equal dtw package")
library(dtw)


test_that("Bub: dtw normalized mv", {
   Q <- matrix(rnorm(100), ncol=2)
   C <- matrix(rnorm(800), ncol=2)
   
   tmp <- IncDTW::dtw(Q, C)
   # tmp$distance
   # tmp$normalized_distance
   
   dist_norm <- tmp$distance/(nrow(Q)+nrow(C))
   dist_norm2 <- tmp$normalized_distance
   expect_equal(dist_norm, dist_norm2)
   
})


test_that("IncDTW:dtw2vec is equal dtw::dtw", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   ws <- NULL
   fooInc <- function(Q,C,ws){ dtw2vec(Q = Q, C = C, ws = ws, step_pattern = "symmetric2") }
   foodtw <- function(Q,C,ws){ dtw::dtw(Q,C, step.pattern = symmetric2,
                                        window.type = "sakoechiba", window.size = ws)}
   
   expect_equal(fooInc(Q, C, ws = NULL)$distance, foodtw(Q, C, ws = Inf)$distance)
   expect_equal(fooInc(Q, C, ws = 50)$normalized_distance,   foodtw(Q, C, ws = 50)$normalizedDistance)
   expect_equal(fooInc(Q, C, ws = 20)$distance,   foodtw(Q, C, ws = 20)$distance)
   expect_equal(fooInc(Q, C, ws = 20)$normalized_distance,   foodtw(Q, C, ws = 20)$normalizedDistance)
})


test_that("IncDTW:dtw2vec is equal dtw::dtw, symmetric2", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   nd <- abs(length(Q) - length(C))# ndiff
   fooVec  <- function(Q,C,ws){ dtw2vec(Q = Q, C = C, ws = ws, step_pattern = "symmetric2")$distance }
   foodtw0 <- function(Q,C,ws){ IncDTW::dtw(Q,C, step_pattern = "symmetric2",
                                        ws = ws)$distance }
   
   foodtw  <- function(Q,C,ws){ dtw::dtw(Q,C, step.pattern = symmetric2,
                                        window.type = "sakoechiba", window.size = ws)$distance }
   
   expect_equal(fooVec(Q, C, ws = NULL), foodtw(Q, C, ws = Inf))
   expect_equal(fooVec(Q, C, ws = nd),   foodtw(Q, C, ws = nd))# fail
   expect_equal(fooVec(Q, C, ws = nd * 3),   foodtw(Q, C, ws = nd * 3))
   
   expect_equal(fooVec(Q, C, ws = NULL), foodtw0(Q, C, ws = NULL))
   expect_equal(fooVec(Q, C, ws = nd),   foodtw0(Q, C, ws = nd))#fail
   expect_equal(fooVec(Q, C, ws = nd * 3),   foodtw0(Q, C, ws = nd * 3))
   
   expect_equal(foodtw(Q, C, ws = Inf),  foodtw0(Q, C, ws = NULL))
   expect_equal(foodtw(Q, C, ws = nd),   foodtw0(Q, C, ws = nd))
   expect_equal(foodtw(Q, C, ws = nd * 3),   foodtw0(Q, C, ws = nd * 3))
   
   if(FALSE){
      fooVec  <- function(Q,C,ws){ dtw2vec(Q = Q, C = C, ws = ws, step_pattern = "symmetric1")$distance }
      foodtw0 <- function(Q,C,ws){ IncDTW::dtw(Q,C, step_pattern = "symmetric1",
                                               ws = ws)$distance }
      
      foodtw  <- function(Q,C,ws){ dtw::dtw(Q,C, step.pattern = dtw::symmetric1,
                                            window.type = "sakoechiba", window.size = ws)$distance }
      counter <- 0
      for(i in 1:100){
         set.seed(i)
         # i <- 3
         Q <- cumsum(rnorm(5))
         C <- Q[3:5] + rnorm(3, 0, 0.5)
         nd <- abs(length(Q) - length(C))# ndiff
         counter <- counter + 
            abs(fooVec(Q, C, ws = nd)- foodtw(Q, C, ws = nd))# fail
            +
               abs(fooVec(Q, C, ws = nd)-   foodtw0(Q, C, ws = nd))            
         
      }
      print(counter)
      
      tmp <- dtw::dtw(Q, C, step.pattern = symmetric2,
                      window.type = "sakoechiba", window.size = nd, keep.internals = T)
      tmp1 <- IncDTW::dtw(Q, C, step_pattern = "symmetric2", ws = nd)
      tmp$costMatrix - tmp1$gcm
      dtw2vec(Q = Q, C = C, ws = nd, step_pattern = "symmetric2")$distance
      
      tmp$costMatrix
      tmp1$gcm
      
   }
   
})



test_that("IncDTW:dtw2vec is equal dtw::dtw, symmetric2_mv", {
   Q <- matrix(rnorm(100), ncol=2)
   C <- matrix(rnorm(80), ncol=2)
   nd <- abs(nrow(Q)-nrow(C))# ndiff
   ws <- 20
   fooVec  <- function(Q,C,ws){ dtw2vec(Q = Q, C = C, ws = ws, step_pattern = "symmetric2", dist_method = "norm2")$distance }
   foodtw0 <- function(Q,C,ws){ dtw2vec_multiv(Q,C, step_pattern = "symmetric2", dist_method = "norm2",
                                            ws = ws)$distance }
   
   foodtw  <- function(Q,C,ws){ dtw::dtw(Q,C, step.pattern = symmetric2,
                                         window.type = "sakoechiba", window.size = ws)$distance }
   
   expect_equal(fooVec(Q, C, ws = NULL), foodtw(Q, C, ws = Inf))
   expect_equal(fooVec(Q, C, ws = 50),   foodtw(Q, C, ws = 50))
   expect_equal(fooVec(Q, C, ws = 20),   foodtw(Q, C, ws = 20))
   
   expect_equal(fooVec(Q, C, ws = NULL), foodtw0(Q, C, ws = NULL))
   expect_equal(fooVec(Q, C, ws = nd),   foodtw0(Q, C, ws = nd))
   expect_equal(fooVec(Q, C, ws = 2*nd), foodtw0(Q, C, ws = 2*nd))
   
   expect_equal(foodtw(Q, C, ws = Inf),  foodtw0(Q, C, ws = NULL))
   expect_equal(foodtw(Q, C, ws = nd),   foodtw0(Q, C, ws = nd))
   expect_equal(foodtw(Q, C, ws = 2*nd), foodtw0(Q, C, ws = 2*nd))
   
   
 
   
})


test_that("Warping paths are equal", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   fooInc <- function(Q,C,ws){ IncDTW::dtw(Q = Q, C = C, ws = ws, step_pattern = 'symmetric1', return_wp = TRUE)$wp }
   foodtw <- function(Q,C,ws){ dtw::dtw(Q,C, step.pattern = symmetric1,
                                        window.type = "sakoechiba", window.size = ws)$stepsTaken }
   
   expect_equal(fooInc(Q, C, ws = NULL), foodtw(Q, C, ws = Inf))
   expect_equal(fooInc(Q, C, ws = 50),   foodtw(Q, C, ws = 50))
   expect_equal(fooInc(Q, C, ws = 20),   foodtw(Q, C, ws = 20))
})



test_that("Equal 1-norm", {
   Q <- cumsum(rnorm(100))
   C <- cumsum(rnorm(100))
   fooInc <- function(Q,C,ws){ IncDTW::dtw2vec(Q = Q, step_pattern = "symmetric1", C = C, ws = ws)$distance }
   foo1norm <- function(Q,C){ sum(abs(Q-C)) }
   
   expect_equal(fooInc(Q, C, ws = 0), foo1norm(Q, C))
})




test_that("Equal dtw2vec and dtw", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   foo2vec <- function(Q,C,ws){ IncDTW::dtw2vec(Q = Q, C = C, ws = ws)$distance }
   foodtw  <- function(Q,C,ws){ IncDTW::dtw(Q, C, ws = ws)$distance }
   foocm <- function(Q, C, ws){ 
      cm <- IncDTW::dtw(Q, C, ws = ws, return_cm = TRUE)$cm
      IncDTW::dtw2vec(Q = cm, C="cm", ws = ws)$distance}
   
   expect_equal(foo2vec(Q, C, ws = NULL), foodtw(Q, C, ws = NULL))
   expect_equal(foo2vec(Q, C, ws = 50),   foodtw(Q, C, ws = 50))
   expect_equal(foo2vec(Q, C, ws = 20),   foodtw(Q, C, ws = 20))
   
   expect_equal(foo2vec(Q, C, ws = NULL), foocm(Q, C, ws = NULL))
   expect_equal(foo2vec(Q, C, ws = 50),   foocm(Q, C, ws = 50))
   expect_equal(foo2vec(Q, C, ws = 20),   foocm(Q, C, ws = 20))
})


test_that("Equal Incremental vs. Scratch vec-based", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   nC <- length(C)
   
   res0  <- idtw2vec(Q = Q, newObs = C[1:(nC-10)], gcm_lc = NULL)#initial
   res1  <- idtw2vec(Q = Q, newObs = C[(nC-9):nC], gcm_lc = res0$gcm_lc_new, nC = nC-10)#incremental
   resfs <- idtw2vec(Q = Q, newObs = C, gcm_lc = NULL)#from scratch
   
   expect_equal(res1$distance, resfs$distance)
   expect_equal(res1$normalized_distance, resfs$normalized_distance)
   
   
   # multivariate
   Q <- matrix(rnorm(20), ncol=2)
   C <- matrix(rnorm(6), ncol=2)
   nC <- nrow(C)
   
   res0  <- idtw2vec(Q = Q, newObs = C[1:(nC-1), ], gcm_lc = NULL)#initial
   res1  <- idtw2vec(Q = Q, newObs = C[(nC-1+1):nC,, drop = FALSE], gcm_lc = res0$gcm_lc_new, nC = nC-1)#incremental
   resfs <- idtw2vec(Q = Q, newObs = C, gcm_lc = NULL)#from scratch
   
   expect_equal(res1$distance, resfs$distance)
   expect_equal(res1$normalized_distance, resfs$normalized_distance)
   
})


test_that("Equal Incremental vs. Scratch vec-based with WS", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   nC <- length(C)
   # skip("skip this: just for testing")
   res0  <- idtw2vec(Q = Q, newObs = C[1:(nC-10)]   , ws = 30, nC = 1    , gcm_lc = NULL)#initial
   res1  <- idtw2vec(Q = Q, newObs = C[(nC-10+1):nC], ws = 30, nC = nC-10, gcm_lc = res0$gcm_lc_new)#incremental
   resfs <-  dtw2vec(Q = Q, C = C,                    ws = 30)#from scratch
   
   expect_equal(res1$distance, resfs$distance)
})



test_that("Equal Incremental vs. Scratch vec-based, with cost matrix", {
   Q <- cumsum(rnorm(100))
   C <- Q[11:100] + rnorm(90, 0, 0.5)
   nC <- length(C)
   # skip("skip this: just for testing")
   cm <- IncDTW::cm(Q, C)
   res00 <- dtw(Q = cm, C = "cm")
   res0  <- idtw2vec_cm(cm = cm[,1:(nC-10), drop = FALSE],     gcm_lc = NULL)#initial
   res1  <- idtw2vec_cm(cm = cm[, (nC-10+1):nC, drop = FALSE], gcm_lc = res0$gcm_lc_new)#incremental
   
   expect_equal(res1$distance, res00$distance)
   
   
   WS <- 30
   res00 <- dtw(Q = cm, C = "cm", ws = WS)
   res0  <- idtw2vec_cm(cm = cm[,1:(nC-10), drop = FALSE], ws = WS,     gcm_lc = NULL)#initial
   res1  <- idtw2vec_cm(cm = cm[, (nC-10+1):nC, drop = FALSE], ws = WS, nC = nC-10, gcm_lc = res0$gcm_lc_new)#incremental
   
   expect_equal(res1$distance, res00$distance)
   expect_equal(res1$normalized_distance, res00$normalized_distance)
   
   
   WS <- 30
   res00 <- dtw(Q = cm, C = "cm", ws = WS)
   res0  <- idtw2vec(Q = cm[,1:(nC-10), drop = FALSE], newObs = "cm", ws = WS,     gcm_lc = NULL)#initial
   res1  <- idtw2vec(Q = cm[, (nC-10+1):nC, drop = FALSE], newObs = "cm",ws = WS, nC = nC-10, gcm_lc = res0$gcm_lc_new)#incremental
   
   expect_equal(res1$distance, res00$distance)
   expect_equal(res1$normalized_distance, res00$normalized_distance)
})






test_that("Equal Incremental vs. Scratch with ws", {
   # skip("test")
   Q <- cumsum(rnorm(10))
   C <- cumsum(rnorm(8))
   WS <- 5
   newObs <-  c(2,3)# new observation
   base <- dtw(Q=Q, C=C, ws = WS, return_diffM = TRUE) # the ordinary calculation
   
   
   #--- recalculation from scratch with new observations
   result0 <- dtw(Q=Q, C=c(C, newObs), ws = WS,  return_diffM = TRUE) # the ordinary calculation
   
   #--- the incremental step with new observations
   result1 <- idtw(Q, C, ws = WS, newObs = newObs, gcm=base$gcm, dm=base$dm, diffM = base$diffM, 
                   return_diffp = TRUE,  return_diffM = TRUE) 
   
   #--- the incremental step with new observations, 
   #     but already calculated additive costMatrix cm_add
   mQ <- matrix(Q, ncol = length(newObs), nrow = length(Q), byrow = FALSE)
   mC <- matrix(newObs, ncol = length(newObs), nrow = length(Q), byrow = TRUE)
   cm_add <- matrix(abs(mQ - mC), ncol = length(newObs))
   result2 <- idtw(Q=cm_add, C="cm_add", ws = WS, newObs = newObs, gcm=base$gcm, dm=base$dm) 
   
   expect_equal(result0$distance, result1$distance)
   expect_equal(result0$distance, result2$distance)
})



test_that("Double Incremental Matrix EQUAL Scratch", {
   # skip("test")
   # in very rare cases the solution of the cheapest path through the gcm is not unique, even though
   # the GCM and DTW distance measure are identical.
   # => the solution of the direction matrix, and the warping path can differ slightly 
   # This happens for example for myseed = 407. 
   # To pass the CRAN tests, I fix the seed.
   myseed <- 123
   set.seed(myseed)
   Q <- cumsum(rnorm(10))
   C <- cumsum(rnorm(8))
   WS <- 5
   newObsQ <-  rnorm(2)# new observation
   newObsC <-  rnorm(3)# new observation
   base <- IncDTW::dtw(Q=Q, C=C, ws = WS, return_diffM = TRUE) # the ordinary calculation
   
   
   #--- recalculation from scratch with new observations
   result0 <- IncDTW::dtw(Q = c(C, newObsC), C = c(Q, newObsQ), ws = WS,  
                          return_diffp = TRUE,return_diffM = TRUE, return_wp = TRUE) # the ordinary calculation
   
   #--- the incremental step with new observations
   result1 <- IncDTW::idtw(Q = Q, C = C, ws = WS, newObs = newObsC, gcm=base$gcm, 
                           dm=base$dm, diffM = base$diffM, 
                           return_diffp = TRUE,  return_diffM = TRUE) 
   
   tmp <- result1$dm
   tmp2 <- which(tmp==2)
   tmp3 <- which(tmp==3)
   tmp[tmp2] <- 3
   tmp[tmp3] <- 2
   tmp <- t(tmp)
   result2 <- IncDTW::idtw(Q = c(C, newObsC), C = Q, ws = WS, newObs = newObsQ, 
                           gcm = t(result1$gcm), dm = tmp, diffM = -t(result1$diffM), 
                           return_diffp = TRUE,  return_diffM = TRUE) 
   
   
   expect_equal(result0$gcm, result2$gcm)
   expect_equal(result0$dm, result2$dm)
   expect_equal(result0$diffp, result2$diffp)
   expect_equal(result0$diffM, result2$diffM)
   expect_equal(result0$wp, result2$wp)
   expect_equal(result0$distance, result2$distance)
})






test_that("Double Incremental Vector EQUAL Scratch", {

   C0 <- cumsum(rnorm(100))
   Q0 <- cumsum(rnorm(80))
   C_new <- cumsum(rnorm(10))
   Q_new <- cumsum(rnorm(10))
   C_update <- c(C0, C_new) 
   Q_update <- c(Q0, Q_new) 
   
   res_scratch <- IncDTW::dtw2vec(Q = Q_update, C = C_update)
   resV_init <- IncDTW::idtw2vec(Q = Q0, newObs = C0, gcm_lc = NULL)
   resV_inc <-  IncDTW::idtw2vec(Q = Q0, newObs = C_new, 
                                gcm_lc = resV_init$gcm_lc_new)
   resV_inc2 <- IncDTW::idtw2vec(Q = C_update, newObs = Q_new, 
                                 gcm_lc = c(resV_init$gcm_lr_new,
                                            resV_inc$gcm_lr_new))
   
   expect_equal(res_scratch$distance, resV_inc2$distance)
   
   
   # again with last row already in first incremntal step
   resV_init <- IncDTW::idtw2vec(Q = Q0, newObs = C0, gcm_lc = NULL)
   resV_inc <-  IncDTW::idtw2vec(Q = Q0, newObs = C_new, 
                                 gcm_lc = resV_init$gcm_lc_new,
                                 gcm_lr = resV_init$gcm_lr_new)
   
   resV_inc2 <- IncDTW::idtw2vec(Q = C_update, newObs = Q_new, 
                                 gcm_lc = resV_inc$gcm_lr_new)
   
   expect_equal(res_scratch$distance, resV_inc2$distance)
   
   
})




test_that("dtw_partial() based on dtw()", {
   
   C <- cos(1:100)
   Q <- c(C, rnorm(100))
   
   tmp <- IncDTW::dtw(Q = Q, C = C, return_QC = T, return_wp = T, step_pattern = "symmetric2")
   ret <- dtw_partial(tmp)
   
   expect_equal(ret$normalized_distance, 0)
   expect_equal(ret$rangeQ[2], 100)
   
   
})




test_that("dtw_partial() based on idtw2vec()", {
   
   C <- cos(1:100)
   Q <- c(C, rnorm(100))
   
   tmp <- IncDTW::idtw2vec(Q = Q, newObs = C, step_pattern = "symmetric2")
   ret <- dtw_partial(tmp)
   
   expect_equal(ret$normalized_distance, 0)
   expect_equal(ret$rangeQ[2], 100)
   
})



test_that("dtw_partial() reverse", {
   
   C <- cos(1:100)
   Q <- c(rnorm(100), C)
   
   tmp <- IncDTW::idtw2vec(Q = rev(Q), newObs = rev(C), step_pattern = "symmetric2")
   #tmp <- dtw(Q = rev(Q), C = rev(C), step_pattern = "symmetric2")
   ret <- dtw_partial(tmp, reverse = TRUE)
   
   expect_equal(ret$normalized_distance, 0)
   expect_equal(ret$rangeQ[1], 101)
   
})


test_that("dec_dm", {
   Q <- cos(1:100)
   C <- cumsum(rnorm(80))
   # the ordinary calculation
   result_base <- dtw(Q=Q, C=C, return_wp = TRUE) 
   
   # the ordinary calculation without the last 4 observations
   result_decr <- dtw(Q=Q, C=C[1:(length(C) - 4)], return_wp = TRUE) 
   # the decremental step: reduce C for 4 observation
   result_decr2 <- dec_dm(result_base$dm, Ndec = 4) 
   
   # compare ii, jj and wp of result_decr and those of 
   expect_equal(result_decr$ii, result_decr2$ii)
   expect_equal(result_decr$jj, result_decr2$jj)
   expect_equal(result_decr$wp, result_decr2$wp)

})