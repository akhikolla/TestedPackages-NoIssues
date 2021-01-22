#'Probabilistic Patient Record Linkage
#'
#'@param data1 either a binary (\code{1} or \code{0} values only) matrix or binary 
#'data frame of dimension \code{n1 x K} whose rownames are the observation identifiers.
#'
#'@param data2 either a binary (\code{1} or \code{0} values only) matrix or a binary
#'data frame of dimension \code{n2 x K} whose rownames are the observation identifiers.
#'
#'@param data1_cont2diff either a matrix or dataframe of continuous features, 
#'such as age, for which the similarity measure uses the difference with 
#'\code{data2_cont2diff}, whose rownames are . Default is \code{NULL}.
#'
#'@param data2_cont2diff either a matrix or dataframe of continuous features, 
#'such as age, for which the similarity measure uses the difference with 
#'\code{data2_cont1diff}, whose rownames are . Default is \code{NULL}.
#'
#'@param eps_plus discrepancy rate between \code{data1} and \code{data2}
#'
#'@param eps_minus discrepancy rate between \code{data2} and \code{data1}
#'
#'@param aggreg_2ways a character string indicating how to merge the posterior two 
#'probability matrices obtained for each of the 2 databases. Four possibility are 
#'currently implemented: \code{"maxnorm"}, \code{"max"}, \code{"min"}, \code{"mean"} 
#'and \code{"prod"}. Default is \code{"mean"}.
#'
#'@param min_prev minimum prevalence for the variables used in matching.
#'Default is 1\%.
#'
#'@param d_max a numeric vector of length \code{K} giving the minimum difference 
#'from which it is considered a discrepancy.
#'
#'@param use_diff logical flag indicating whether continuous differentiable variables should be used in the 
#'
#'@param dates1 matrix or dataframe of dimension \code{n1 x K} including the concatenated dates intervals for each corresponding 
#'diagnosis codes in \code{data1}. Default is \code{NULL} in which case dates are not used.
#'
#'@param dates2 matrix or dataframe of dimension \code{n2 x K} including the concatenated dates intervals for each corresponding 
#'diagnosis codes in \code{data2}. Default is \code{NULL} in which case dates are not used. See details.
#'
#'@details \code{Dates:} the use of \code{dates1} and \code{dates2} requires that at least one date interval matches across 
#'\code{dates1} and \code{dates2} for claiming an agreement on a diagnosis code between \code{data1} and \code{data2}, 
#'in addition of having that very same code recorded in both.
#'
#@param eps_inf discrepancy rate for the differentiable variables
#'
#'@references Hejblum BP, Weber G, Liao KP, Palmer N, Churchill S, Szolovits P, Murphy S, Kohane I, Cai T 
#'Probabilistic Record Linkage of De-Identified Research Datasets Using Diagnosis Codes, \emph{submitted}, 2017.
#'
#'@importFrom landpred VTM
#'@importFrom fGarch dsstd sstdFit
#'
#'@return a matrix of size \code{n1 x n2} with the posterior probability of matching for each \code{n1*n2} pair
#'
#'@examples
#'set.seed(123)
#'ncodes <- 500
#'npat <- 200
#'incid <- abs(rnorm(n=ncodes, 0.15, 0.07))
#'bin_codes <- rbinom(n=npat*ncodes, size=1,  prob=rep(incid, npat))
#'bin_codes_mat <- matrix(bin_codes, ncol=ncodes, byrow = TRUE)
#'data1_ex <- bin_codes_mat[1:(npat/2+npat/10),]
#'data2_ex <- bin_codes_mat[c(1:(npat/10), (npat/2+npat/10 + 1):npat), ]
#'rownames(data1_ex) <- paste0("ID", 1:(npat/2+npat/10), "_data1")
#'rownames(data2_ex) <- paste0("ID", c(1:(npat/10), (npat/2+npat/10 + 1):npat), "_data2")
#'
#'if(interactive()){
#'res <- recordLink(data1 = data1_ex, data2 = data2_ex, 
#'                  use_diff = FALSE, eps_minus = 0.01, eps_plus = 0.01)
#'round(res[c(1:3, 19:23), c(1:3, 19:23)], 3)
#'}
#'
#'@export

recordLink <- function(data1, data2, dates1 = NULL, dates2 = NULL,
                       eps_plus, eps_minus, aggreg_2ways = "mean", 
                       min_prev = 0.01, 
                       data1_cont2diff = NULL, data2_cont2diff = NULL,
                       #eps_inf, 
                       d_max, use_diff=TRUE){
  
  datesNotNull <- (!is.null(dates1) & !is.null(dates2))
  nb_feat <- ncol(data1)
  if(ncol(data2)!=nb_feat){stop("Number of columns in data2 is different from data1")}
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  if((is.null(data1_cont2diff) | is.null(data2_cont2diff)) && use_diff){
    stop("cannot 'use_diff' when 'data1_cont2diff' and/or 'data2_cont2diff' is NULL\n Probably need to set 'use_diff = FALSE'")
  }
  if((is.null(dates1) & !is.null(dates2)) | (!is.null(dates1) & is.null(dates2))){
    stop("missing one of the dates tables")
  }
  
  if(use_diff){
    ndiff <- ncol(data1_cont2diff) 
    stopifnot(ncol(data2_cont2diff)==ndiff)
    
    #if(length(eps_inf)==1 && ndiff>1){
    #  eps_inf <- rep(eps_inf, ndiff)
    #}
    if(length(d_max)==1 && ndiff>1){
      d_max <- rep(d_max, ndiff)
    }
    #stopifnot(length(eps_inf)==ndiff)
    stopifnot(length(d_max)==ndiff)
  }
  
  ind1 <- rownames(data1)
  ind2 <- rownames(data2)
  
  freq_select <- (colMeans(data1) > min_prev)
  data1_bin <- data1[,freq_select, drop=FALSE]
  data2_bin <- data2[,freq_select, drop=FALSE]
  if(datesNotNull){
    dates1_fs <- dates1[,freq_select, drop=FALSE]
    dates2_fs <- dates2[,freq_select, drop=FALSE]
  }
  #rm(list=c("data1", "data2"))
  #n_freq10 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.10 | colSums(data1_bin)/nrow(data1_bin)>0.90)
  #n_freq5 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.05 | colSums(data1_bin)/nrow(data1_bin)>0.95)
  #n_freq1 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.01 | colSums(data1_bin)/nrow(data1_bin)>0.99)
  
  #removing info too rare in any of the 2 database ?
  col_in <- (pmax(apply(X=(data1_bin==0), MARGIN=2, FUN=mean), apply(X=data2_bin==0, MARGIN=2, FUN=mean)) < 0.995 &  
               pmin(apply(data1_bin==0,2,mean),apply(data2_bin==0,2,mean)) > 0.005)
  if(sum(col_in)<1){stop("No features whith adequate prevalence")}
  
  data1_bin <- data1_bin[, col_in, drop=FALSE]
  data2_bin <- data2_bin[, col_in, drop=FALSE]
  if(datesNotNull){
    dates1_fs_ci <- dates1_fs[, col_in, drop=FALSE]
    dates2_fs_ci <- dates2_fs[, col_in, drop=FALSE]
  }
  
  pi1 <- apply(data1_bin==1, MARGIN=2, FUN=mean)
  pi2 <- apply(data2_bin==1, MARGIN=2, FUN=mean)
  K0 <- ncol(data1_bin)
  if(length(eps_plus)==1){
    eps_p <- rep(eps_plus, K0)
  }else if(length(eps_plus)==nb_feat){
    eps_p <- eps_plus[freq_select][col_in]
  }else{
    stop("Length of 'eps_plus' doesn't make sense")
  }
  
  if(length(eps_minus)==1){
    eps_n <- rep(eps_minus, K0)
  }else if(length(eps_minus)==nb_feat){
    eps_n <- eps_minus[freq_select][col_in]
  }else{
    stop("Length of 'eps_minus' doesn't make sense")
  }
  
  
  # need matrices for the C function
  data1_bin <- as.matrix(data1_bin)
  data2_bin <- as.matrix(data2_bin)
  if(datesNotNull){
    dates1_fs_ci <- as.matrix(dates1_fs_ci)
    dates2_fs_ci <- as.matrix(dates2_fs_ci)
  }
  zeros_p <- which(eps_p==0)
  if(length(zeros_p)>0){
    eps_p[zeros_p] <- eps_p[zeros_p] + 1E-6
  }
  zeros_n <- which(eps_n==0)
  if(length(zeros_n)>0){
    eps_n[zeros_n] <- eps_n[zeros_n] + 1E-6
  }
  keep_eps <- which(eps_p!=1 & eps_n!=1)

  if(datesNotNull){
    dist_bin <- loglikC_bin_wDates(Bmat = data2_bin[, keep_eps], Amat = data1_bin[, keep_eps],
                                   Bdates = dates2_fs_ci[, keep_eps], Adates = dates1_fs_ci[, keep_eps],
                                   eps_p = eps_p[keep_eps], eps_n = eps_n[keep_eps], 
                                   piA = pi1[keep_eps], piB = pi2[keep_eps])
  }else{
    dist_bin <- loglikC_bin(Bmat = data2_bin[, keep_eps], Amat = data1_bin[, keep_eps], eps_p = eps_p[keep_eps], 
                             eps_n = eps_n[keep_eps], piA = pi1[keep_eps], piB = pi2[keep_eps])
  }
  
  if(use_diff){
    pi_inf <- numeric(ndiff)
    for(j in 1:ndiff){
      pi_inf[j] <- mean((unlist(lapply(data1_cont2diff[, j, drop=FALSE], function(x){x-data2_cont2diff[, j, drop=FALSE]}))<d_max[j]))
    }
    
    dist_diff <- loglikratioC_diff_arbitrary(Bmat = as.matrix(t(data2_cont2diff)), 
                                             Amat =  as.matrix(t(data1_cont2diff)), 
                                             d_max = d_max, 
                                             cost = rep(-1E+4, ndiff)
    )
    
    
    dist_all <- dist_bin + dist_diff
  }else{
    dist_all <- dist_bin
  }
  rm(list=c("data1_bin", "data2_bin"))
  colnames(dist_all) <- ind2
  rownames(dist_all) <- ind1
  
  #sstdFit generates warnings - to be ignored
  if(length(dist_all)>10000){
    sstdFit_1way <-  try(suppressWarnings(fGarch::sstdFit(x = sample(dist_all, 10000))), silent=TRUE)
    if(inherits(sstdFit_1way, "try-error")){sstdFit_1way <-  try(suppressWarnings(fGarch::sstdFit(x = sample(dist_all,10000))), silent=TRUE)}
    tmp_est <-  sstdFit_1way$estimate
    for(i in 1:4){
      sstdFit_1way_sub <-  try(suppressWarnings(fGarch::sstdFit(x = sample(dist_all, 10000))), silent=TRUE)
      if(inherits(sstdFit_1way_sub, "try-error")){sstdFit_1way_sub <-  try(suppressWarnings(fGarch::sstdFit(x = sample(dist_all, 10000))), silent=TRUE)}
      tmp_est <-  rbind(tmp_est, sstdFit_1way_sub$estimate)
    }
    sstdFit_1way$est <- colMeans(tmp_est, na.rm = TRUE)
  }else{
    sstdFit_1way <-  try(suppressWarnings(fGarch::sstdFit(x = as.vector(dist_all))), silent=TRUE)
    if(inherits(sstdFit_1way, "try-error")){sstdFit_1way <-  try(suppressWarnings(fGarch::sstdFit(x = as.vector(dist_all))), silent=TRUE)}
    if(inherits(sstdFit_1way, "try-error")){
      stop("Error in fitting the Student-t distribution")
    }
    sstdFit_1way$est <- sstdFit_1way$estimate
  }
  
  inter_x <- 0.1
  obs_dist_x_1way <- seq(min(dist_all), max(dist_all), by=inter_x)
  
  
  ### select cutoff
  # # empirical solution
  # # single skew-t parametric estimation
  # fitskewt_dens_est <- fGarch::dsstd(obs_dist_x_1way,
  #                                    mean=sstdFit_1way$est[1], 
  #                                    sd=sstdFit_1way$est[2],
  #                                    nu=sstdFit_1way$est[3],
  #                                    xi=sstdFit_1way$est[4])
  # peak_fitskewt_dens_est <- max(fitskewt_dens_est)
  # rho1_hat_empirical_1way <- min((obs_dist_x_1way[obs_dist_x_1way >= peak_fitskewt_dens_est])[fitskewt_dens_est[obs_dist_x_1way >= peak_fitskewt_dens_est]/peak_fitskewt_dens_est < 1/min(n1,n2)])
  
  # analytical solution
  nu <-  sstdFit_1way$est[3]
  s <-  sstdFit_1way$est[2]
  m <-  sstdFit_1way$est[1]
  xi <-  sstdFit_1way$est[4]
  Beta <- exp(lgamma(0.5) - lgamma((nu+1)/2) + lgamma(nu/2))
  m1 <- 2*sqrt(nu-2)/((nu-1)*Beta)
  mu <- m1*(xi-1/xi)
  sigma <- sqrt((1-m1^2)*(xi^2+1/xi^2)+2*m1^2-1)
  d1max_x <- m-(mu-xi*sqrt((nu-2)/(nu+2)))*s/sigma  
  
  d1 <- function(x){
    #c <- 2*sigma*gamma((nu+1)/2)/(s*(xi+1/xi)*sqrt(pi*nu-2)*gamma(nu/2))
    c <- exp(log(2) + log(sigma) + lgamma((nu+1)/2) - log(s*(xi+1/xi)*sqrt(pi*nu-2)) - lgamma(nu/2))
    d <- -c*(nu+1)*sigma/(xi^2*(nu-2)*s)
    y <- (x-m)/s*sigma+mu
    return(d*y*(1+y^2/(xi^2*(nu-2)))^(-(nu+3)/2))
  }
  
  d2 <- function(x){
    #c <- 2*sigma*gamma((nu+1)/2)/(s*(xi+1/xi)*sqrt(pi*nu-2)*gamma(nu/2))
    c <- exp(log(2) + log(sigma) + lgamma((nu+1)/2) - log(s*(xi+1/xi)*sqrt(pi*nu-2)) - lgamma(nu/2))
    d <- -c*(nu+1)*sigma/(xi^2*(nu-2)*s)
    y <- (x-m)/s*sigma+mu
    return(d*sigma/s*(1+y^2/(xi^2*(nu-2)))^(-(nu+5)/2)*(1-y^2*(nu+2)/xi^2*(nu-2)))
  }
  
  rho1_hat <- obs_dist_x_1way[obs_dist_x_1way>d1max_x][min(which(abs(d1(obs_dist_x_1way[obs_dist_x_1way>d1max_x]))<(1/min(n1,n2)) & abs(d2(obs_dist_x_1way[obs_dist_x_1way>d1max_x]))<(1/min(n1,n2))))]
  if(is.na(rho1_hat)){
    rho1_hat <- obs_dist_x_1way[obs_dist_x_1way>d1max_x][min(which(abs(d1(obs_dist_x_1way[obs_dist_x_1way>d1max_x]))<(1/min(n1,n2))))]
  }
  if(is.na(rho1_hat)){
    rho1_hat <- 0
  }
  rm(list="obs_dist_x_1way")
  
  
  dist_all_select <- dist_all > rho1_hat # rho1_hat_empirical_1way
  nmatch_hat_1way <- min(sum(colSums(dist_all_select)>=1), sum(rowSums(dist_all_select)>=1)) # Empirical Bayes prior on P(M=1)
  rm(list="dist_all_select")
  
  prop_match_1way <- nmatch_hat_1way/(n1*n2)
  
  rank_match_1 <- matchProbs_rank_full_C(dist_all, prop_match=prop_match_1way)
  rownames(rank_match_1) <- ind1
  colnames(rank_match_1) <- ind2
  rank_match_2 <- t(matchProbs_rank_full_C(t(dist_all), prop_match=prop_match_1way))
  rm(list="dist_all")
  rownames(rank_match_2) <- ind1
  colnames(rank_match_2) <- ind2
  mode(rank_match_1) <- "numeric" 
  mode(rank_match_2) <- "numeric"
  
  if(aggreg_2ways == "maxnorm"){
    rank_match_add0_1 <- cbind(rank_match_1, 1-rowSums(rank_match_1))
    probMax_1 <- apply(rank_match_add0_1, MARGIN=1, FUN=max)
    rm(list="rank_match_add0_1")
    rank_match_add0_2 <- rbind(rank_match_2, 1-colSums(rank_match_2))
    probMax_2 <- apply(rank_match_add0_2, MARGIN=2, FUN=max)
    rm(list="rank_match_add0_2")
  }
  
  probMatch_mat <- switch(aggreg_2ways,
                          maxnorm = rank_match_1*rank_match_2/matrix(probMax_1, ncol=1)%*%matrix(probMax_2, nrow=1),
                          min = pmin(rank_match_1, rank_match_2),
                          max = pmax(rank_match_1, rank_match_2),
                          mean = (rank_match_1+rank_match_2)/2,
                          prod = rank_match_1*rank_match_2
  )
  
  
  return(probMatch_mat)
  
}


