pair_fast_mgm = function(data,         # n x p data matrix
                type,         # p vector indicating the type of each variable
                level,        # p vector indivating the levels of each variable
                lambdaSeq,    # sequence of considered lambda values (default to glmnet default)
                lambdaSel,    # way of selecting lambda: CV vs. EBIC
                lambdaFolds,  # number of folds if lambdaSel = 'CV'
                lambdaGam,    # EBIC hyperparameter gamma, if lambdaSel = 'EBIC'
                alphaSeq,     # sequence of considered alpha values (elastic net), default = 1 = lasso
                alphaSel,     # way of selecting lambda: CV vs. EBIC
                alphaFolds,   # number of folds if alphaSel = 'CV'
                alphaGam,     # EBIC hyperparameter gamma, if alphaSel = 'EBIC',
                
                ruleReg,      # rule to combine d+1 neighborhood estimates (defaults to 'AND')
                weights,      # p vector of observation weights for weighted regression
                threshold,    # defaults to 'LW', see helpfile
                method,       # glm vs. 'linear'; for now only implement glm
                binarySign,   # should a sign be computed for binary models (defaults to NO)
                scale,        # should gaussian variables be scaled? defaults to TRYE
                verbatim,     # turns off all notifications
                pbar,         #
                warnings,     #
                overparameterize, # if TRUE, uses the over-parameterized version,
                thresholdCat,
                signInfo,
                ...
){


###### pairwise(k=2) fast mgm #####
k=2
p <- ncol(data)
n <- nrow(data)

colnames(data)[1:p] <- paste("V", 1:p, ".", sep = "")

# ----- Fill in Defaults -----

if(missing(lambdaSeq)) lambdaSeq <- NULL
if(missing(lambdaSel)) lambdaSel <- 'CV'
if(missing(lambdaFolds)) lambdaFolds <- 10
if(missing(lambdaGam)) lambdaGam <- .25
if(missing(alphaSeq)) alphaSeq <- 1
if(missing(alphaSel)) alphaSel <- 'CV'
if(missing(alphaFolds)) alphaFolds <- 10
if(missing(alphaGam)) alphaGam <- .25
if(missing(ruleReg)) ruleReg <- 'AND'
if(missing(weights)) weights <- rep(1, n)
if(missing(threshold)) threshold <- 'LW'
if(missing(method)) method <- 'glm'
if(missing(binarySign)) binarySign <- FALSE
if(missing(scale)) scale <- TRUE
if(missing(verbatim)) verbatim <- FALSE
if(missing(pbar)) pbar <- TRUE
if(missing(warnings)) warnings <- TRUE
if(missing(overparameterize)) overparameterize <- FALSE
if(missing(signInfo)) signInfo <- TRUE
if (missing(thresholdCat)) 
  if (overparameterize) 
    thresholdCat <- TRUE
else thresholdCat <- TRUE

if(verbatim) pbar <- FALSE
if(verbatim) warnings <- FALSE

# Switch all warnings off
if(!warnings) {
  oldw <- getOption("warn")
  options(warn = -1)
}


d <- k - 1
emp_lev <- rep(NA, p)
ind_cat <- which(type == "c")
if (length(ind_cat) > 0) 
  for (i in 1:length(ind_cat)) emp_lev[ind_cat][i] <- length(unique(data[, 
                                                                         ind_cat[i]]))
emp_lev[which(type != "c")] <- 1

if(!missing(level)) {
  # Check whether provided levels are equal to levels in the data
  level_check <- level != emp_lev
  if(sum(level_check) > 0) stop(paste0('Provided levels not equal to levels in data for variables ',paste((1:p)[level_check], collapse = ', ')))
  # if not provided, do nothing, because the argument is not actually necessary
}
# Normalize weights (necessary to ensure that nadj makes sense)
if(!missing(weights)) weights <- weights / max(weights)
nadj <- sum(weights) # calc adjusted n



level <- emp_lev
nadj <- sum(weights)
ind_Gauss <- which(type == "g")



# Scale Gaussians
ind_Gauss <- which(type == 'g')
if(scale) for(i in ind_Gauss) data[, i] <- scale(data[, i])


# ----- Basic Checks -----

if(!(threshold %in% c('none', 'LW', 'HW'))) stop('Please select one of the three threshold options "HW", "LW" and "none" ')

if(nrow(data) < 2) ('The data matrix has to have at least 2 rows.')

if(missing(data)) stop('No data provided.')
if(k<2) stop('The order of interactions should be at least k = 2 (pairwise interactions)')
if(ncol(data)<3) stop('At least 3 variables required')
if(missing(type)) stop('No type vector provided.')

if(sum(!(type %in% c('g', 'c', 'p')))>0) stop("Only Gaussian 'g', Poisson 'p' or categorical 'c' variables permitted.")
if(any(is.na(data))) stop('No missing values permitted.')
if(any(!is.finite(as.matrix(data)))) stop('No infinite values permitted.')


if(any(!(apply(data, 2, class) %in% c('numeric', 'integer')))) stop('Only integer and numeric values permitted.')

if(ncol(data) != length(type)) stop('Number of variables is not equal to length of type vector.')
if(!missing(level)) if(ncol(data) != length(level)) stop('Number of variables is not equal to length of level vector.')
if(nrow(data) != length(weights)) stop('Number of observations is not equal to length of weights vector.')

# Are Poisson variables integers?
if('p' %in% type) {
  ind_Pois <- which(type == 'p')
  nPois <- length(ind_Pois)
  v_PoisCheck <- rep(NA, length=nPois)
  for(i in 1:nPois) v_PoisCheck[i] <- sum(data[, ind_Pois[i]] != round(data[, ind_Pois[i]])) > 0
  if(sum(v_PoisCheck) > 0) stop('Only integers permitted for Poisson variables.')
}

# ----- Checking glmnet minimum Variance requirements -----

glmnetRequirements(data = data, type = type, weights = weights)
ind_cat <- which(type == "c")
ind_binary <- rep(NA, length(ind_cat))
ind_binary <- as.logical(ind_binary)
if (length(ind_cat) > 0) {
  for (i in 1:length(ind_cat)) ind_binary[i] <- length(unique(data[, 
                                                                   ind_cat[i]])) == 2
}
if (sum(ind_binary) > 0) {
  check_binary <- rep(NA, sum(ind_binary))
  for (i in 1:sum(ind_binary)) check_binary[i] <- sum(!(unique(data[, 
                                                                    ind_cat[ind_binary][i]]) %in% c(0, 1)))
  if (binarySign) {
    if (sum(check_binary) > 0) 
      stop(paste0("If binarySign = TRUE, all binary variables have to be coded {0,1}. Not satisfied in variable(s) ", 
                  paste(ind_cat[ind_binary][check_binary > 0], 
                        collapse = ", ")))
  }
}
mgmobj <- list(call = NULL, pairwise = list(wadj = NULL, 
                                            signs = NULL, edgecolor = NULL), factorgraph = list(graph = NULL, 
                                                                                                signs = NULL, edgecolor = NULL, order = NULL), rawfactor = list(indicator = NULL, 
                                                                                                                                                                weightsAgg = NULL, weights = NULL, signs = NULL), intercepts = NULL, 
               nodemodels = list())
mgmobj$call <- list(data = NULL, type = type, level = level, 
                    lambdaSeq = lambdaSeq, lambdaSel = lambdaSel, lambdaFolds = lambdaFolds, 
                    lambdaGam = lambdaGam, alphaSeq = alphaSeq, alphaSel = alphaSel, 
                    alphaFolds = alphaFolds, alphaGam = alphaGam, k = k, 
                    ruleReg = ruleReg, weights = weights, threshold = threshold, 
                    method = method, binarySign = binarySign, scale = scale, 
                    verbatim = verbatim, pbar = pbar, warnings = warnings, 
                     overparameterize = overparameterize)
data <- as.data.frame(data)
for (i in which(type == "c")) data[, i] <- as.factor(data[, 
                                                          i])
if (pbar == TRUE) 
  pb <- txtProgressBar(min = 0, max = p, initial = 0, 
                       char = "-", style = 3)
npar_standard <- rep(NA, p)



for (v in 1:p) {
  if (d > (p - 1)) {
    stop("Order of interactions cannot be larger than the number of predictors.")
  }
  else if (d == 1) {
    form <- as.formula(paste(colnames(data)[v], "~ (.)"))
  }
  else {
    form <- as.formula(paste(colnames(data)[v], "~ (.)^", 
                             d))
  }
  X_standard <- model.matrix(form, data = data)[, -1]
  npar_standard[v] <- ncol(X_standard)
  if (overparameterize) {
    X_over <- ModelMatrix(data = data[, -v], type = type[-v], 
                          level = level[-v], labels = colnames(data)[-v], 
                          d = d)
    X <- X_over
  }
  else {
    X <- X_standard
  }
  y <- as.numeric(data[, v])
  n_alpha <- length(alphaSeq)
  if (alphaSel == "CV") {
    l_alphaModels <- list()
    ind <- sample(1:alphaFolds, size = n, replace = TRUE)
    v_mean_OOS_deviance <- rep(NA, n_alpha)
    if (n_alpha > 1) {
      for (a in 1:n_alpha) {
        l_foldmodels <- list()
        v_OOS_deviance <- rep(NA, alphaFolds)
        for (fold in 1:alphaFolds) {
          train_X <- X[ind != fold, ]
          train_y <- y[ind != fold]
          test_X <- X[ind == fold, ]
          test_y <- y[ind == fold]
          n_train <- nrow(train_X)
          nadj_train <- sum(weights[ind != fold])
          l_foldmodels[[fold]] <- nodeEst(y = train_y, 
                                          X = train_X, lambdaSeq = lambdaSeq, lambdaSel = lambdaSel, 
                                          lambdaFolds = lambdaFolds, lambdaGam = lambdaGam, 
                                          alpha = alphaSeq[a], weights = weights[ind != 
                                                                                   fold], n = n_train, nadj = nadj_train, 
                                          v = v, type = type, level = level, emp_lev = emp_lev, 
                                          overparameterize = overparameterize, thresholdCat=thresholdCat)
          LL_model <- calcLL(X = test_X, y = test_y, 
                             fit = l_foldmodels[[fold]]$fitobj, type = type, 
                             level = level, v = v, weights = weights[ind == 
                                                                       fold], lambda = l_foldmodels[[fold]]$lambda, 
                             LLtype = "model")
          LL_saturated <- calcLL(X = test_X, y = test_y, 
                                 fit = l_foldmodels[[fold]]$fitobj, type = type, 
                                 level = level, v = v, weights = weights[ind == 
                                                                           fold], lambda = l_foldmodels[[fold]]$lambda, 
                                 LLtype = "saturated")
          v_OOS_deviance[fold] <- 2 * (LL_saturated - 
                                         LL_model)
        }
        v_mean_OOS_deviance[a] <- mean(v_OOS_deviance)
      }
      alpha_select <- alphaSeq[which.min(v_mean_OOS_deviance)]
    }
    else {
      alpha_select <- alphaSeq
    }
    model <- nodeEst(y = y, X = X, lambdaSeq = lambdaSeq, 
                     lambdaSel = lambdaSel, lambdaFolds = lambdaFolds, 
                     lambdaGam = lambdaGam, alpha = alpha_select, 
                     weights = weights, n = n, nadj = nadj, v = v, 
                     type = type, level = level, emp_lev = emp_lev, 
                     overparameterize = overparameterize,thresholdCat=thresholdCat)
    mgmobj$nodemodels[[v]] <- model
  }
  if (alphaSel == "EBIC") {
    l_alphaModels <- list()
    EBIC_Seq <- rep(NA, n_alpha)
    for (a in 1:n_alpha) {
      l_alphaModels[[a]] <- nodeEst(y = y, X = X, 
                                    lambdaSeq = lambdaSeq, lambdaSel = lambdaSel, 
                                    lambdaFolds = lambdaFolds, lambdaGam = lambdaGam, 
                                    alpha = alphaSeq[a], weights = weights, n = n, 
                                    nadj = nadj, v = v, type = type, level = level, 
                                    emp_lev = emp_lev, overparameterize = overparameterize,thresholdCat=thresholdCat)
      EBIC_Seq[a] <- l_alphaModels[[a]]$EBIC
    }
    ind_minEBIC_model <- which.min(EBIC_Seq)
    mgmobj$nodemodels[[v]] <- l_alphaModels[[ind_minEBIC_model]]
  }
  if (pbar == TRUE) 
    setTxtProgressBar(pb, v)
}



########## begin to modify #######

CalculateTau=function(mat,npar_standard,threshold){
  mat=mat[-1,]
  if (threshold == "LW") 
    tau <- sqrt(d) * sqrt(sum(mat^2)) * 
      sqrt(log(npar_standard)/n)
  if (threshold == "HW") 
    tau <- d * sqrt(log(npar_standard)/n)
  if (threshold == "none") 
    tau <- 0
  mat[abs(mat) < tau] = 0
  return (as.matrix(mat) )
}

CalculateCoef=function(pair1,pair2){
  t=level[-pair1]-1
  t[t==0]=1
  npar_standard=sum(t)
  model_obj=mgmobj$nodemodels[[pair1]]$model
  
  if (type[pair2]=='c'){
    a=paste0(nodename[pair2],1:(level[pair2]-1) )
    if (is.list(model_obj)) {
      
      model_obj_ni=sapply(1:length(model_obj), function(x) CalculateTau(model_obj[[x]],npar_standard,threshold)   )
      rownames(model_obj_ni)=rownames(model_obj[[1]])[-1]
      coef=sapply(1:ncol(model_obj_ni), function(x) model_obj_ni[a,x]   )
      
    } 
    else {
      
      model_obj=CalculateTau(model_obj,npar_standard,threshold)
      coef=model_obj[a,1]
    }
  } else {
    a=nodename[pair2]
    if (is.list(model_obj)) {
      model_obj_ni=sapply(1:length(model_obj), function(x) CalculateTau(model_obj[[x]],npar_standard,threshold)   )
      rownames(model_obj_ni)=rownames(model_obj[[1]])[-1]
      coef=sapply(1:ncol(model_obj_ni), function(x) model_obj_ni[a,x]   )
    }
    else {
      model_obj=CalculateTau(model_obj,npar_standard,threshold)
      coef=model_obj[a,1]
    }
  }
  return(mean( abs(coef) ))
}

alledge=which( upper.tri(matrix(1,nrow=p,ncol=p))!=0,arr.ind = T)
nodename=names(data)

wadj=matrix(0,p,p)
for (i in 1:nrow(alledge)) {
  coef1=CalculateCoef(alledge[i,1],alledge[i,2])
  coef2=CalculateCoef(alledge[i,2],alledge[i,1])
  paircoef=c(coef1,coef2)
  
  if (ruleReg == "AND") 
    parcompr = mean( paircoef ) * (!(0 %in% paircoef))
  if (ruleReg == "OR") 
    parcompr = mean( paircoef )
  wadj[alledge[i,1],alledge[i,2]]=parcompr
}

wadj=wadj+t(wadj)

return(wadj)

}




environment(pair_fast_mgm)=asNamespace('mgm')



