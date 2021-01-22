ypinterim.default <- function(time, event, group, spendfun, critvalue = NULL, repnum = 1E4,
                                  bound = 50, seed.fix = 0, ...) {
  
  free <- ifelse(is.logical(seed.fix) & (seed.fix == FALSE), TRUE, FALSE)

  if((free == FALSE) & is.numeric(seed.fix)) {
    if (exists(".Random.seed", .GlobalEnv)) {
      oldseed <- .Random.seed
    } else {
      oldseed <- NULL
    }
    set.seed(seed.fix)
  }
  
  nmonitoring <- length(spendfun)
  len_usercrt <- length(critvalue)
  display  <- (len_usercrt + 1):nmonitoring
  alpha_spen <- diff(c(0, spendfun))
  
  if(!exists("time")) stop("'time' is missing.")
  if(!exists("event")) stop("'event' is missing.")
  if(!exists("group")) stop("'group' is missing.")
  if(!exists("spendfun")) stop("'spendfun' is missing.")
  
  if(!all(spendfun == cummax(spendfun))) stop("The elements of 'spendfun' must be increasing.")
  if(any(nmonitoring != c(ncol(time), ncol(event)))) stop("The number of looks of 'time', 'event' and 'spendfun' do not match.")
  if(len_usercrt >= nmonitoring) stop("'critvalue' should be at least one less than the number of the looks.")
  
  y_iv <- as.matrix(time)
  n <- nrow(y_iv)
  
  for (i_mon in 1:nmonitoring) {
    each_Y_iv <- y_iv[, i_mon]
    if (length(unique(each_Y_iv)) != length(each_Y_iv)) {
      y_iv[, i_mon] <- each_Y_iv + runif(n, 0, 1e-50)
    }
  }
  
  delta_iv <- as.matrix(event) == 1
  z_i <- as.matrix(group) == 1
  
  o_all <- apply(y_iv, 2, function(x) {
    order(x)
  })
  
  quan1 <- quan2 <- quan_s <- v1 <- v2 <- v_s <- matrix(NA, nrow = n, ncol = nmonitoring)
  fitbeta_all <- NULL
  
  for (o_i in 1:nmonitoring) {
    oy_iv <- y_iv[o_all[, o_i], o_i]
    od_iv <- delta_iv[o_all[, o_i], o_i]
    oz_i <- z_i[o_all[, o_i]]
    fitbeta <- fun_estimate(oy = oy_iv, od = od_iv, oz = oz_i, bound = bound)
    fitbeta_all <- cbind(fitbeta_all, fitbeta)
    gamma1 <- exp(-oz_i * fitbeta[1])
    gamma2 <- exp(-oz_i * fitbeta[2])
    K <- n:1
    Lambda1 <- cumsum(od_iv * gamma1/K)
    Lambda2 <- cumsum(od_iv * gamma2/K)
    P <- exp(-Lambda2)
    PL <- c(1, P[1:(n - 1)])
    R <- cumsum(PL * od_iv * gamma1/K)/P
    
    tt1 <- exp(fitbeta[1])
    tt2 <- exp(fitbeta[2])
    
    psi1 <- (1 + R)/(1/tt1 + 1/tt2 * R)
    psi2 <- 1/psi1
    
    yl <- K
    yzl <- cumsum(oz_i[n:1])[n:1]
    
    quan1[, o_i] <- od_iv * psi1 * (oz_i - yzl/yl)
    quan2[, o_i] <- od_iv * psi2 * (oz_i - yzl/yl)
    
    v1[, o_i] <- psi1^2 * (yzl/yl - (yzl/yl)^2) * od_iv
    v2[, o_i] <- psi2^2 * (yzl/yl - (yzl/yl)^2) * od_iv
    
    quan_s[, o_i] <- od_iv * (oz_i - yzl/yl)
    v_s[, o_i] <- (yzl/yl - (yzl/yl)^2) * od_iv
  }
  
  sd1_den <- sqrt(colSums(v1))
  sd2_den <- sqrt(colSums(v2))
  if (nmonitoring > 1) sds_den <- apply(event, 2, sum)

  w1 <- colSums(quan1)/sd1_den
  w2 <- colSums(quan2)/sd2_den

  testvalue <- apply(cbind(w1, w2), 1, FUN = function(x) max(abs(x)))
  
  c_names <- paste("c", 1:nmonitoring, sep = "")
  pcrt_names <- paste("pcrt", 1:nmonitoring, sep = "")
  pstat_names <- paste("pstat", 1:nmonitoring, sep = "")
  
  if (!is.null(critvalue)) {
    for (cc in 1:nmonitoring) {
      if (cc <= len_usercrt) {
        assign(c_names[cc], critvalue[cc])
        assign(pcrt_names[cc], NA)
        assign(pstat_names[cc], NA)
      } else {
        assign(c_names[cc], NULL)
        assign(pcrt_names[cc], NULL)
        assign(pstat_names[cc], NULL)
      }
    }
  } else {
    for (cc in 1:nmonitoring) {
      assign(c_names[cc], NULL)
      assign(pcrt_names[cc], NULL)
      assign(pstat_names[cc], NULL)
    } 
  }
  
  for (c_i in 1:nmonitoring) {
    
      if(alpha_spen[c_i] < 1E-4) {
        
      prob_lowertail <- function(crt) {
        if (c_i == 1) res <- log(2*pnorm(crt, lower.tail = FALSE))
        if (c_i == 2) {
          rho <- sqrt(sds_den[1]/sds_den[2])
          joint_f <- function(z) {
            dnorm(z)*(pnorm((crt - z*rho)/sqrt(1 - rho^2)) - pnorm((-crt - z*rho)/sqrt(1 - rho^2)))
          }
          joint_a <- integrate(f = joint_f, lower = -get(c_names[1]), upper = get(c_names[1]))
          res <- log(1 - joint_a$value)
        }
        if (c_i > 2) {
          sig <- outer(sds_den[1:c_i], sds_den[1:c_i], FUN = function(x,y) sqrt(x/y))
          sig[lower.tri(sig)] <- sig[upper.tri(sig)]
          simz <- MASS::mvrnorm(n = 1000000, mu = rep(0, c_i), Sigma = sig)
          satisfied <- apply(simz[, -c_i], 1, function(x) all(abs(x) <= unlist(mget(c_names[1:(c_i-1)], inherits = TRUE))))
          res <- log(mean(abs(simz[satisfied, c_i]) > crt))
        }
        return(res)
      }
      
      if (is.null(get(c_names[c_i]))) {
        crt <- fun_interpol(alphas = alpha_spen[c_i], f = prob_lowertail)
      } else {
        crt <- get(c_names[c_i])
      }
      
      pcrt <- 2*pnorm(crt, lower.tail = FALSE)  
      
      if (testvalue[c_i] > qnorm(1 - 1E-4/2)) {
        pstat <- 2*pnorm(testvalue[c_i], lower.tail = FALSE)  
      } else {
        mquan1 <- as.matrix(quan1[, 1:c_i])
        mquan2 <- as.matrix(quan2[, 1:c_i])
        t_tilde <- awlrstat(repnum = repnum, n = n, quan1 = mquan1,
                            quan2 = mquan2, sd1 = sd1_den[1:c_i], sd2 = sd2_den[1:c_i])
        pstat <- mean(t_tilde[,c_i] - testvalue[c_i] > 0)
      }
    } else {
      mquan1 <- as.matrix(quan1[, 1:c_i])
      mquan2 <- as.matrix(quan2[, 1:c_i])
      
      t_tilde <- awlrstat(repnum = repnum, n = n, quan1 = mquan1,
                          quan2 = mquan2, sd1 = sd1_den[1:c_i], sd2 = sd2_den[1:c_i])
      
      if (c_i == 1) {
        target_quantile <- t_tilde[, 1]
      } else {
        ind <- 1
        for (con_i in 1:(c_i - 1)) {
          ind <- (t_tilde[, con_i] <= get(c_names[con_i])) * ind
        }
        target_quantile <- t_tilde[as.logical(ind), c_i]
      }
      
      if (is.null(get(c_names[c_i]))) {
        sort_target <- sort(target_quantile)
        crt <- sort_target[length(sort_target) - floor(repnum*alpha_spen[c_i])]
      } else {
        crt <- get(c_names[c_i])
      }

      pcrt <- mean(t_tilde[, c_i] - crt > 0)
      
      if (testvalue[c_i] > qnorm(1 - 1E-4/2)) {
        pstat <- 2*pnorm(testvalue[c_i], lower.tail = FALSE)  
      } else {
        pstat <- mean(t_tilde[,c_i] > testvalue[c_i])
      }
    }
    
    assign(c_names[c_i], crt)
    assign(pcrt_names[c_i], pcrt)
    assign(pstat_names[c_i], pstat)
    
  }
  
  c_values <- unlist(mget(c_names))
  pcrt_values <- unlist(mget(pcrt_names))
  pstat_values <- unlist(mget(pstat_names))

  colnames(fitbeta_all) <- 1:nmonitoring
  rownames(fitbeta_all) <- c("b1", "b2")
  
  
  if (!is.null(oldseed)) {
    on.exit( { .Random.seed <<- oldseed} )
  } else {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  
  result <- list()
  result$teststat <- testvalue
  result$pvalue <- pstat_values
  result$all.critvalue <- c_values
  result$new.critvalue <- c_values[display]
  result$user.critvalue <- critvalue
  result$pvalue.bnd <- pcrt_values
  result$beta <- fitbeta_all
  result$bound <- bound
  result$repnum <- repnum
  result$spendfun <- spendfun
  result$dspendfun <- diff(c(0, spendfun))
  result$n.looks <- nmonitoring 
  result$call <- match.call()

  class(result) <- "ypinterim"
  
  result
}

