#' Association testing by combining several matching thresholds
#' 
#' Computes association test p-values from a generalized linear model for each considered 
#' threshold, and computes a p-value for the combination of all the envisioned thresholds 
#' through Fisher's method using perturbation resampling.
#'
#'@param match_prob matching probabilities matrix (e.g. obtained through \code{\link{recordLink}}) of 
#'dimensions \code{n1 x n2}.
#'
#'@param y response variable of length \code{n1}. Only binary phenotypes are supported at the moment.
#'
#'@param x a \code{matrix} or a \code{data.frame} of predictors of dimensions \code{n2 x p}. 
#'An intercept is automatically within the function.
#'
#'@param thresholds a vector (possibly of length \code{1}) containing the different threshold 
#'to use to call a match. Default is \code{seq(from = 0.5, to = 0.95, by = 0.05)}.
#'
#'@param nb_perturb the number of perturbation used for the p-value combination.
#'Default is 200.
#'
#'@param dist_family a character string indicating the distribution family for the glm. 
#'Currently, only \code{'gaussian'} and  \code{'binomial'} are supported. Default 
#'is \code{'gaussian'}.
#'
#'@param impute_strategy a character string indicating which strategy to use to impute x 
#'from the matching probabilities \code{match_prob}. Either \code{"best"} (in which 
#'case the highest probable match above the threshold is imputed) or \code{"weighted average"}
#'(in which case weighted mean is imputed for each individual who has at least
#'one match with a posterior probability above the threshold). Default is 
#'\code{"weighted average"}.
#'
#'@importFrom landpred VTM
#'@importFrom fGarch dsstd sstdFit
#'@importFrom stats binomial glm na.omit rnorm as.formula model.matrix
#'
#'@return a list containing the following:
#'\itemize{
#'   \item \code{influencefn_pvals} p-values obtained from influence function perturbations
#'   with the covariates as columns and the \code{thresholds} as rows, with an additional row 
#'   at the top for the combination 
#'   \item \code{wald_pvals} a matrix containing the p-values obtained from the Wald 
#'   test with the covariates as columns and the \code{thresholds} as rows
#'   \item \code{ptbed_pvals} a list containing, for each covariates, a matrix with
#'   the \code{nb_perturb} perturbed p-values with the different \code{thresholds}
#'   as rows
#'   \item \code{theta_impute} a matrix of the estimated coefficients from the glm when imputing 
#'   the weighted average for covariates (as columns) with the \code{thresholds} as rows
#'   \item \code{sd_theta} a matrix of the estimated SD (from the influence function) of the 
#'   coefficients from the glm when imputing the weighted average for covariates (as columns),
#'   with the \code{thresholds} as rows
#'   \item \code{ptbed_theta_impute} a list containing, for each covariates, a matrix with
#'   the \code{nb_perturb} perturbed estimated coefficients from the glm when imputing 
#'   the weighted average for covariates, with the different \code{thresholds}
#'   as rows
#'   \item \code{impute_strategy} a character string indicating which impute 
#'   strategy was used (either \code{"weighted average"} or \code{"best"})
#'}
#'
#'@export
#'
#'@examples
#'#rm(list=ls())
#'res <- list()
#'n_sims <- 1#5000
#'for(n in 1:n_sims){
#'x <- matrix(ncol=2, nrow=99, stats::rnorm(n=99*2))
#'
#'#plot(density(rbeta(n=1000, 1,2)))
#'match_prob <- matrix(rbeta(n=103*99, 1, 2), nrow=103, ncol=99)
#'
#'
#'#y <- rnorm(n=103, 1, 0.5)
#'#res[[n]] <- test_combine(match_prob, y, x, dist_family="gaussian")$influencefn_pvals
#'y <- rbinom(n=103, 1, prob=0.5)
#'res[[n]] <- test_combine(match_prob, y, x, dist_family="binomial")$influencefn_pvals
#'cat(n, "/", n_sims, "\n", sep="")
#'}
#'size <- matrix(NA, ncol=nrow(res[[1]]), nrow=ncol(res[[1]])-2)
#'colnames(size) <- rownames(res[[1]])
#'rownames(size) <- colnames(res[[1]])[-(-1:0 + ncol(res[[1]]))]
#'for(i in 1:(ncol(res[[1]])-2)){
#'  size[i, ] <- rowMeans(sapply(res, function(m){m[, i]<0.05}), na.rm = TRUE)
#'}
#'size
#'

test_combine <- function(match_prob, y, x,
                     thresholds = seq(from = 0.5, to = 0.95, by = 0.05), #c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
                     nb_perturb = 200,
                     dist_family = c("gaussian", "binomial"),
                     impute_strategy = c("weighted average", "best")){
  
  if(length(dist_family)>1){
    dist_family <- dist_family[1]
  }
  
  # sanity checks
  stopifnot(is.matrix(match_prob))
  
  if(length(which(is.na(x))) > 0){
    warning("x contains NA/nan: to be able to provide results associated observations will be removed")
    x_toremove <- unique(which(is.na(x), arr.ind = TRUE)[, "row"])
    x <- x[-x_toremove, ]
    match_prob <- match_prob[, -x_toremove]
  }
  if(is.data.frame(x)){
    x <- stats::model.matrix(stats::as.formula(paste0("~", paste(colnames(x), collapse=" + "))), data = x)[, -1, drop = FALSE]
  }else if(is.matrix(x)){
    if(all(x[,1]==1) | ifelse(!is.null(colnames(x)), colnames(x)[1] == "(Intercept)", FALSE)){
      x <- x[, -1, drop = FALSE]
    }
  }else{
    stop("x is neither a data.frame nor a matrix")
  }
  stopifnot(is.vector(y))
  
  
  n1 <- length(y)
  n2 <- nrow(x)
  stopifnot(nrow(match_prob) == n1)
  stopifnot(ncol(match_prob) == n2)
  
  if(length(which(is.na(y))) > 0){
    warning("y contains NA/nan: to be able to provide results associated observations will be removed")
    y_toremove <- which(is.na(y))
    y <- y[-y_toremove]
    match_prob <- match_prob[-y_toremove, ]
  }
  n1 <- length(y)
  stopifnot(nrow(match_prob) == n1)
  
  nb_thresholds <- length(thresholds)
  mismatch_avg <- numeric(nb_thresholds)
  names(mismatch_avg) <- thresholds
  
  if(!(dist_family %in% c("gaussian", "binomial"))){
    stop("'gaussian' or 'binomial' are the only valid values for dist_family currently supported")
  }
  
  
  stopifnot(is.character(impute_strategy))
  if(length(impute_strategy) > 1){
    impute_strategy <- impute_strategy[1]
  }
  stopifnot(impute_strategy %in% c("weighted average", "best"))
  
  # initializing results
  eta <- list()
  ncoef <- ncol(x) +1
  theta_avg <- matrix(NA, ncol = ncoef, nrow = length(thresholds))
  rownames(theta_avg) <- as.character(thresholds)
  wald_pvals <- matrix(NA, ncol = ncoef, nrow = length(thresholds))
  rownames(wald_pvals) <- as.character(thresholds)
  
  for(i in 1:length(thresholds)){
    cut_p <- thresholds[i]
    prob_sup_cut <- (match_prob > cut_p)
    
    #construct the data-frame for the glm
    xi <- rowSums(prob_sup_cut) > 0
    n_rho <- sum(xi)
    
    y_match <- y*xi
    
    match_prob_sel <-  diag(1*xi) %*% match_prob
    if(impute_strategy == "best"){
      xi_NA <- xi
      xi_NA[!xi] <- NA
      x_impute <- diag(1*xi_NA) %*% x[max.col(match_prob_sel), ]
    }else if(impute_strategy == "weighted average"){
      x_impute <- match_prob_sel %*% x / rowSums(match_prob_sel) #diag(1/rowSums(match_prob_sel)) %*% match_prob_sel %*% x/rowSums(match_prob_sel)
    }else{
      stop("'strategy' is neither 'best' nor 'weighted average'")
    }
    
    #y_match_sub <- y[xi]
    #x_best_sub <- x[max.col(match_prob[xi, ]), ]
    #x_impute_sub <- match_prob[xi, ]%*%x/rowSums(match_prob[xi, ])
    impute_fit_summary <- summary(stats::glm(y_match ~ x_impute, family = dist_family, na.action = stats::na.omit))
    theta_avg[i, ] <- impute_fit_summary$coef[, "Estimate", drop=FALSE]
    if(dist_family == "binomial"){
      wald_pvals[i, ] <- impute_fit_summary$coef[, "Pr(>|z|)", drop=FALSE]
    }else if(dist_family == "gaussian"){
      wald_pvals[i, ] <- impute_fit_summary$coef[, "Pr(>|t|)", drop=FALSE]
    }else{
      stop("dist_family is neither 'gaussian' nor 'binomial'")
    }
    # Z_sub <- stats::model.matrix( ~ x_impute_sub)
    # I_rho <- 1/n_rho*crossprod(apply(Z_sub, 2, function(colu){colu*sqrt(expit_dev1(Z_sub%*%theta))}))
    # eta[[as.character(cut_p)]] <- 1/n_rho*solve(I_rho)%*%t(Z_sub)%*%diag(x=(y_match_sub - expit(Z_sub%*%theta)[, "Estimate"]))
    # sqrt(apply(eta[[as.character(cut_p)]], 1, crossprod))
    
    x_impute_noNA <-  x_impute 
    x_impute_noNA[is.na(x_impute[, 1])] <- 0
    Z <- diag(xi) %*% stats::model.matrix( ~ x_impute_noNA)
    
    if(dist_family == "binomial"){
      I_rho <- 1/n_rho*crossprod(apply(Z, 2, function(colu){colu*sqrt(expit_dev1(Z %*% theta_avg[i, ]))}))
      eta[[as.character(cut_p)]] <- 1/n_rho*solve(I_rho) %*% t(Z) %*% diag(x=(y_match - xi*expit(Z %*% theta_avg[i, ])[, 1]))
    }else if(dist_family == "gaussian"){
      I_rho <- 1/n_rho*crossprod(Z)
      eta[[as.character(cut_p)]] <- 1/n_rho*solve(I_rho) %*% t(Z) %*% diag(x=(y_match - xi*(Z %*% theta_avg[i, ])[, 1]))
    }else{
      stop("dist_family is neither 'gaussian' nor 'binomial'")
    }
    
    # fit <- stats::glm(y_match ~ x_impute, family = stats::binomial, na.action = stats::na.omit)
    # summary(fit)$coef
    # sqrt(apply(eta[[as.character(cut_p)]], 1, crossprod))
    # solve(vcov(fit))/n_rho
  }
  colnames(theta_avg) <- rownames(impute_fit_summary$coefficients)
  colnames(wald_pvals) <- rownames(impute_fit_summary$coefficients)
  
  sigma <- list()
  sds <- matrix(ncol=ncol(theta_avg), nrow = length(thresholds))
  colnames(sds) <- colnames(theta_avg)
  rownames(sds) <- thresholds
  eta_mat <- list()
  for(j in 1:ncol(theta_avg)){
    eta_mat[[j]] <- sapply(X = eta, FUN = function(m){m[j, ]})
    sigma[[colnames(theta_avg)[j]]] <- crossprod(eta_mat[[j]])  
    sds[, j] <- sqrt(diag(sigma[[colnames(theta_avg)[j]]]))
  }
  names(eta_mat) <- colnames(theta_avg)
  
  
  # Combine p-values
  pvals <- pval_zscore(theta_avg, sds)
  fisher_comb <- apply(pvals, MARGIN = 2, FUN = comb_pvals)
  
  B <- nb_perturb #number of perturbations
  theta_avg_star <- lapply(eta_mat, function(m){crossprod(m, matrix(stats::rnorm(n1*B, mean = 0, sd = 1), nrow = n1, ncol = B))})
  pvals_star <- list()
  for(k in 1:ncol(sds)){
    pvals_star[[k]] <- apply(theta_avg_star[[k]], MARGIN = 2, FUN = pval_zscore, sigma = sds[, k])
  }
  if(length(thresholds)>1){
    fisher_comb_star <- lapply(pvals_star, function(m){apply(m, MARGIN = 2, FUN = comb_pvals)})
  }else{
    fisher_comb_star <- lapply(pvals_star, FUN = function(v){sapply(v, comb_pvals)})
  }
  pval_combined <- matrix(nrow=1, ncol=ncol(pvals))
  row.names(pval_combined) <- "Combined p-value"
  colnames(pval_combined) <- colnames(pvals)
  for(k in 1:ncol(sds)){
    pval_combined[1, k] <- 1-sum(fisher_comb[k] >= fisher_comb_star[[k]])/B
  }

  return(list("influencefn_pvals" = cbind.data.frame(rbind(pval_combined, pvals), 
                                                     "matching_threshold" = c("combined", rownames(pvals)),
                                                     "impute_strategy" = impute_strategy),
              "wald_pvals" = wald_pvals,
              "ptbed_pvals" = pvals_star,
              "theta_impute" = theta_avg,
              "sd_theta"=sds,
              "ptbed_theta_impute" = theta_avg_star,
              "impute_strategy" = impute_strategy)
  )
  
}


