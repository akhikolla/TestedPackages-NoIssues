#' Create Function-based Historical Matrix Projection Model
#' 
#' Function \code{flefko3()} returns function-based historical MPMs corresponding
#' to the patches and years given, including the associated component transition
#' and fecundity matrices, data frames detailing the characteristics of the
#' ahistorical stages used and historical stage pairs created, and a data frame
#' characterizing the patch and year combinations corresponding to these matrices.
#' Unlike \code{\link{rlefko3}()}, this function currently does not currently
#' distinguish populations.
#'
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{all} if
#' matrices should be estimated for all patches. Defaults to \code{all}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for each
#' potentially new individual (row) born to each reproductive stage (column).
#' Entries act as multipliers on fecundity, with 1 equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function
#' describing transitions to be overwritten either with given values or 
#' with other estimated transitions.
#' @param data The original historical demographic data frame used to 
#' estimate vital rates (class \code{hfvdata}). The original data frame is 
#' required in order to initialize years and patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate models.
#' If given, then \code{surv_model}, \code{obs_model}, \code{size_model}, \code{repst_model}, \code{fec_model},
#' \code{jsurv_model}, \code{jobs_model}, \code{jsize_model}, \code{jrepst_model}, \code{paramnames},
#' \code{yearcol}, and \code{patchcol} are not required. One or more of these models should
#' include size or reproductive status in time \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can 
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any survival 
#' probability model given in \code{modelsuite}. This model must have been developed 
#' in a modeling exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then 
#' will overwrite any observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing the impacts 
#' of times \emph{t} and \emph{t}-1.
#' @param size_model A linear model predicting size. This can be a model of class
#' \code{glm} or \code{glmer}, both of which require a predicted poisson variable under a 
#' log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian response is 
#' assumed. If given, then will overwrite any size model given in \code{modelsuite}.  
#' This model must have been developed in a modeling exercise testing the impacts 
#' of times \emph{t} and \emph{t}-1.
#' @param repst_model A linear model predicting reproduction probability. This 
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any reproduction 
#' probability model given in \code{modelsuite}.  This model must have been developed 
#' in a modeling exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param fec_model A linear model predicting fecundity. This can be a model
#' of class \code{glm} or \code{glmer}, and requires a predicted poisson variable under a 
#' log link. If given, then will overwrite any fecundity model given in 
#' \code{modelsuite}. This model must have been developed in a modeling exercise 
#' testing the impacts of times \emph{t} and \emph{t}-1.
#' @param jsurv_model A linear model predicting juvenile survival probability. This
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted 
#' binomial variable under a logit link. If given, then will overwrite any juvenile
#' survival probability model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any juvenile observation probability model given in \code{modelsuite}. This model
#' must have been developed in a modeling exercise testing the impacts of times \emph{t}
#' and \emph{t}-1.
#' @param jsize_model A linear model predicting juvenile size. This can be a model
#' of class \code{glm} or \code{glmer}, both of which require a predicted poisson variable 
#' under a log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian 
#' response is assumed. If given, then will overwrite any juvenile size model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in the previous year. This can be a model 
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial variable under a 
#' logit link. If given, then will overwrite any reproduction probability model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param paramnames A dataframe with two columns, the first showing the general
#' model terms that will be used in matrix creation, and the second showing the
#' equivalent terms used in modeling. Only required if \code{modelsuite} is not 
#' supplied.
#' @param inda A numeric value to use for individual covariate a. Defaults to 0.
#' @param indb A numeric value to use for individual covariate b. Defaults to 0.
#' @param indc A numeric value to use for individual covariate c. Defaults to 0.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile reproduction probability.
#' @param repmod A scalar multiplier of fecundity. Defaults to 1.
#' @param yearcol The variable name or column number corresponding to year 
#' in time \emph{t} in the dataset. Not needed if \code{modelsuite} is supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing time step coefficients 
#' are set to 0.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing patch coefficients 
#' are set to 0.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing time step and patch coefficients, if either \code{year.as.random} or
#' \code{patch.as.random} is set to TRUE. Defaults to 0.
#' @param negfec A logical value denoting whether fecundity values estimated to be
#' negative should be reset to 0. Defaults to FALSE.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated solely with 0 transitions. These are only removed in cases where the
#' associated row and column sums in ALL matrices estimated equal 0. Defaults to
#' FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return an
#' object of class \code{lefkoMat}, which is a list with the following elements:
#' 
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' years. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame showing the patch and year of each matrix in order.
#' In \code{flefko3()}, only one population may be analyzed at once, and so \code{pop = NA}.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite} input.}
#'
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", "Sz5nr",
#'                  "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", "Sz4r",
#'                  "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                          obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                          indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                            individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                            sizeacol = "lnVol88", repstracol = "Intactseed88",
#'                            fecacol = "Intactseed88", deadacol = "Dead1988",
#'                            nonobsacol = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathrepmln <- matrix(0, 21, 21)
#' lathrepmln[1, c(13:21)] <- 0.345
#' lathrepmln[2, c(13:21)] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"),
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' lathmodelsln3 <- modelsearch(lathvertln, historical = TRUE, approach = "mixed", suite = "main", 
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl",bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2",year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE, quiet = TRUE)
#' 
#' lathmat3ln <- flefko3(year = "all", patch = "all", stageframe = lathframeln, 
#'                       modelsuite = lathmodelsln3, data = lathvertln, 
#'                       repmatrix = lathrepmln, overwrite = lathover3,
#'                       patchcol = "patchid", yearcol = "year2", year.as.random = FALSE,
#'                       patch.as.random = FALSE, reduce = FALSE)
#'                       
#' summary(lathmat3ln)
#' }
#' 
#' @export
flefko3 <- function(year = "all", patch = "all", stageframe, repmatrix = NA, overwrite = NA,
                    data = NA, modelsuite = NA, surv_model = NA, obs_model = NA, size_model = NA, 
                    repst_model = NA, fec_model = NA, jsurv_model = NA, jobs_model = NA, 
                    jsize_model = NA, jrepst_model = NA, paramnames = NA, inda = 0, indb = 0,
                    indc = 0, surv_dev = 0, obs_dev = 0, size_dev = 0, repst_dev = 0, fec_dev = 0, 
                    jsurv_dev = 0, jobs_dev = 0, jsize_dev = 0, jrepst_dev = 0, repmod = 1, 
                    yearcol = NA, patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE, 
                    randomseed = 0, negfec = FALSE, reduce = FALSE) {
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or equivalents supplied either through the modelsuite option or through the paramnames input parameter.")
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {stop("Need original vertical dataset to set proper limits on year and patch.", call. = FALSE)}
  if (!any(class(data) == "data.frame")) {stop("Need original vertical dataset used in modeling to proceed.", call. = FALSE)}
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element(tolower(year), "all")) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.", call. = FALSE)
  }
  
  if (is.na(patch) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch option when using a modelsuite in which patch is designated.")
  }
  
  if (is.character(patchcol)) {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element(tolower(patch), "all")) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (all(is.na(repmatrix))) {
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (!is.numeric(inda) | !is.numeric(indb) | !is.numeric(indc)) {
    stop("Individual covariate values must equal valid numbers.", call. = FALSE)
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$bin_size_ctr)))))) {
    stop("Function flefko3() requires size to be numeric rather than categorical.", call. = FALSE)
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$stage_id <- as.numeric(stageframe$stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  stageframe$fullstage <- apply(as.matrix(c(1:dim(stageframe)[1])), 1, function(X) {
    paste(stageframe$bin_size_ctr[X], stageframe$repstatus[X])
  })
  
  instages <- length(stageframe$stage_id)
  
  overwrite <- .overwrite_reassess(overwrite, stageframe, historical = TRUE)
  
  # Next the data frame carrying all raw values and element indices for matrix element estimation
  allstages.list <- theoldpizzle(stageframe, overwrite, repmatrix, finalage = 0, style = 0, cont = 0)
  allstages <- do.call("cbind.data.frame", allstages.list)
  
  maxsize <- max(c(allstages$a.size3, allstages$a.size2n, allstages$a.size2o, allstages$a.size1), na.rm = TRUE)
  
  allstages <- allstages[(which(allstages$b.index321 != -1)),]
  
  # Now we work up the models
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
  }
  
  set.seed(randomseed)
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      rvarssummed <- sum(size_proxy$variances[,"vcov"])
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else {
    sizedist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  # This creates a list of pop, patch, and year in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame(NA, X, as.matrix(year));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- NA
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    
    listofyears <- cbind.data.frame(NA, NA, as.matrix(year))
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- NA
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  # A few extra tidbits required for the core matrix estimator to work
  total.matrix.dim <- (length(stageframe$bin_size_ctr) - 1)^2
  
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  # The next call runs the core matrix estimator function and creates all matrices
  madsexmadrigal <- lapply(yearlist, jerzeibalowski, allstages, surv_proxy, obs_proxy, size_proxy, repst_proxy,
                           fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jrepst_proxy, inda, indb, indc, 
                           surv_dev, obs_dev, size_dev, repst_dev, fec_dev, jsurv_dev, jobs_dev, jsize_dev, 
                           jrepst_dev, total.matrix.dim, repmod, rvarssummed, sigma, jrvarssummed, jsigma, 
                           maxsize, sizedist, fecdist, negfec)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  pairings1 <- expand.grid(stage_id_2 = stageframe$stage_id[1:(dim(stageframe)[1] - 1)], 
                           stage_id_1 = stageframe$stage_id[1:(dim(stageframe)[1] - 1)])
  pairings2 <- expand.grid(stage_2 = stageframe$stage[1:(dim(stageframe)[1] - 1)], 
                           stage_1 = stageframe$stage[1:(dim(stageframe)[1] - 1)])
  hstages <- cbind.data.frame(pairings1, pairings2)
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  if (reduce == TRUE) {
    drops <- .reducer3(a_list, u_list, f_list, hstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    hstages <- drops$hstages
  }
  
  rownames(hstages) <- c(1:dim(hstages)[1])
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages, ahstages = ahstages, 
                 labels = listofyears[,c(1:3)], matrixqc = qcoutput1, modelqc = qcoutput2)
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Reduce Matrix Dimensions By Eliminating Empty Stages
#' 
#' \code{.reducer3()} identifies empty stages in a set of historical matrices and
#' removes them from all matrices. It is used within \code{\link{flefko3}()} and 
#' \code{\link{rlefko3}()}.
#' 
#' @param A List of population projection matrices, from a \code{lefkoMat} object.
#' @param U List of surviva-transition matrices corresponding to \code{A}.
#' @param F List of fecundity matrices corresponding to \code{A}.
#' @param hstages Data frame giving the names and identities of historical stage 
#' pairs used to create matrices.
#' 
#' @return Returns a list of reduced \code{A}, \code{U}, and \code{F} matrices, plus the reduced
#' \code{hstages} object.
#' 
#' @keywords internal
#' @noRd
.reducer3 <- function(A, U, F, hstages) {
  stagepatterns <- lapply(A, function(X) {
    matrix.sums <- colSums(X) + rowSums(X)
    return(matrix.sums)
  })
  
  used.stages.mat <- do.call("rbind", stagepatterns)
  used.stages.ovr <- colSums(used.stages.mat)
  keep.stages <- which(used.stages.ovr > 0)
  
  Ared <- lapply(A, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Ured <- lapply(U, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Fred <- lapply(F, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  hstred <- hstages[keep.stages,]
  
  return(list(A = Ared, U = Ured, F = Fred, hstages = hstred))
}

#' Create Function-based Ahistorical Matrix Projection Model
#'
#' Function \code{flefko2()} returns ahistorical MPMs corresponding to the patches
#' and years given, including the associated component transition and fecundity
#' matrices, a data frame detailing the characteristics of the ahistorical stages
#' used, and a data frame characterizing the patch and year combinations 
#' corresponding to these matrices. Unlike \code{\link{rlefko2}()} and \code{\link{rlefko3}()},
#' this function currently does not currently distinguish populations.
#'
#' @param year A variable corresponding to year or observation time, or a set of
#' such values, given in values associated with the year term used in linear model
#' development. Can also equal \code{all}, in which case matrices will be estimated
#' for all years. Defaults to \code{all}.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{all} if
#' matrices should be estimated for all patches. Defaults to \code{all}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for each
#' potentially new individual (row) born to each reproductive stage (column).
#' Entries act as multipliers on fecundity, with 1 equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function
#' describing transitions to be overwritten either with given values or 
#' with other estimated transitions.
#' @param data The original historical demographic data frame used to estimate
#' vital rates (class \code{hfvdata}). The original data frame is required in order
#' to initialize years and patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate models.
#' If given, then \code{surv_model}, \code{obs_model}, \code{size_model}, \code{repst_model}, \code{fec_model},
#' \code{jsurv_model}, \code{jobs_model}, \code{jsize_model}, \code{jrepst_model}, \code{paramnames},
#' \code{yearcol}, and \code{patchcol} are not required. No models should include size or
#' reproductive status in time \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can 
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any survival 
#' probability model given in \code{modelsuite}. This model must have been developed 
#' in a modeling exercise testing only the impacts of time \emph{t}.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any observation probability model given in \code{modelsuite}. This model must
#' have been developed in a modeling exercise testing only the impacts of time \emph{t}.
#' @param size_model A linear model predicting size. This can be a model of class
#' \code{glm} or \code{glmer}, both of which require a predicted poisson variable under a 
#' log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian response is 
#' assumed. If given, then will overwrite any size model given in \code{modelsuite}.  
#' This model must have been developed in a modeling exercise testing only the
#' impacts of time \emph{t}.
#' @param repst_model A linear model predicting reproduction probability. This can
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any reproduction 
#' probability model given in \code{modelsuite}.  This model must have been developed 
#' in a modeling exercise testing only the impacts of time \emph{t}.
#' @param fec_model A linear model predicting fecundity. This can be a model of
#' class \code{glm} or \code{glmer}, and requires a predicted poisson variable under a 
#' log link. If given, then will overwrite any fecundity model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing only the
#' impacts of time \emph{t}.
#' @param jsurv_model A linear model predicting juvenile survival probability. This
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any 
#' juvenile survival probability model given in \code{modelsuite}. This model must 
#' have been developed in a modeling exercise testing only the impacts of time \emph{t}.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any juvenile observation probability model given in \code{modelsuite}. This model
#' must have been developed in a modeling exercise testing only the 
#' impacts of time \emph{t}.
#' @param jsize_model A linear model predicting juvenile size. This can be a model
#' of class \code{glm} or \code{glmer}, both of which require a predicted poisson variable 
#' under a log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian 
#' response is assumed. If given, then will overwrite any juvenile size model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing only the impacts of time \emph{t}.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in the previous year. This can be a model 
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial variable under a 
#' logit link. If given, then will overwrite any reproduction probability model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing only the impacts of time \emph{t}.
#' @param paramnames A dataframe with two columns, the first showing the general
#' model terms that will be used in matrix creation, and the second showing the
#' equivalent terms used in modeling. Only required if \code{modelsuite} is not 
#' supplied.
#' @param inda A numeric value to use for individual covariate a. Defaults to 0.
#' @param indb A numeric value to use for individual covariate b. Defaults to 0.
#' @param indc A numeric value to use for individual covariate c. Defaults to 0.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile reproduction probability.
#' @param repmod A scalar multiplier of fecundity. Defaults to 1.
#' @param yearcol The variable name or column number corresponding to year in time
#' \emph{t} in the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param patchcol The variable name or column number corresponding to patch in the
#' dataset. Not needed if a \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing time step coefficients are
#' set to 0.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing patch coefficients are
#' set to 0.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing time step and patch coefficients, if either \code{year.as.random} or
#' \code{patch.as.random} is set to TRUE. Defaults to 0.
#' @param negfec A logical value denoting whether fecundity values estimated to be
#' negative should be reset to 0. Defaults to FALSE.
#' @param reduce A logical value denoting whether to remove ahistorical stages
#' associated solely with 0 transitions. These are only removed in cases where the
#' associated row and column sums in ALL matrices estimated equal 0. Defaults to
#' FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return
#' either an object of class \code{lefkoMat}. Output includes:
#'
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' years. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{Null for ahistorical matrices.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the patch and year of each matrix in order.
#' In \code{flefko2()}, only one population may be analyzed at once, and so \code{pop = NA}}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the modelsuite input.}
#' 
#' Please note that this function will yield incorrect estimates if the models
#' utilized incorporate state in time \emph{t}-1. Only use models developed testing
#' ahistorical effects.
#'
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", "Sz5nr",
#'                  "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", "Sz4r",
#'                  "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                          obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                          indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                            individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                            sizeacol = "lnVol88", repstracol = "Intactseed88",
#'                            fecacol = "Intactseed88", deadacol = "Dead1988",
#'                            nonobsacol = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathrepmln <- matrix(0, 21, 21)
#' lathrepmln[1, c(13:21)] <- 0.345
#' lathrepmln[2, c(13:21)] <- 0.054
#' 
#' lathover2 <- overwrite(stage3 = c("Sd", "Sdl"), stage2 = c("Sd", "Sd"),
#'                        givenrate = c(0.345, 0.054))
#' 
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "mixed", suite = "main",
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl", bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE, quiet = TRUE)
#'                              
#' lathmat2ln <- flefko2(year = "all", patch = "all", stageframe = lathframeln, 
#'                       modelsuite = lathmodelsln2, data = lathvertln, 
#'                       repmatrix = lathrepmln, overwrite = lathover2,
#'                       patchcol = "patchid", yearcol = "year2",
#'                       year.as.random = FALSE, patch.as.random = FALSE, 
#'                       reduce = FALSE)
#' 
#' summary(lathmat2ln)
#' }
#' 
#' @export
flefko2 <- function(year = "all", patch = "all", stageframe, repmatrix = NA, overwrite = NA,
                    data = NA, modelsuite = NA, surv_model = NA, obs_model = NA, size_model = NA, 
                    repst_model = NA, fec_model = NA, jsurv_model = NA, jobs_model = NA, 
                    jsize_model = NA, jrepst_model = NA, paramnames = NA, inda = 0, indb = 0,
                    indc = 0, surv_dev = 0, obs_dev = 0, size_dev = 0, repst_dev = 0, fec_dev = 0, 
                    jsurv_dev = 0, jobs_dev = 0, jsize_dev = 0, jrepst_dev = 0, repmod = 1, yearcol = NA, 
                    patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE, 
                    randomseed = 0, negfec = FALSE, reduce = FALSE) {
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or equivalents supplied either through the modelsuite option or through the paramnames input parameter.")
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {stop("Need original vertical dataset to set proper limits on year and patch.", call. = FALSE)}
  if (!any(class(data) == "data.frame")) {stop("Need original vertical dataset used in modeling to proceed.", call. = FALSE)}
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element(tolower(year), "all")) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.", call. = FALSE)
  }
  
  if (is.na(patch) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch option when using a modelsuite in which patch is designated.")
  }
  
  if (is.character(patchcol)) {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element(tolower(patch), "all")) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (all(is.na(repmatrix))) {    #This bit needs to be redone to check for fecundity values in the dataset and use them to determine the proper reproductive stages
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (!is.numeric(inda) | !is.numeric(indb) | !is.numeric(indc)) {
    stop("Individual covariate values must equal valid numbers.", call. = FALSE)
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$size)))))) {
    stop("Function flefko2() requires size to be numeric rather than categorical.")
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$stage_id <- as.numeric(stageframe$stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  stageframe$fullstage <- apply(as.matrix(c(1:dim(stageframe)[1])), 1, function(X) {
    paste(stageframe$bin_size_ctr[X], stageframe$repstatus[X])
  })
  
  instages <- length(stageframe$stage_id)
  
  overwrite <- .overwrite_reassess(overwrite, stageframe, historical = FALSE)
  
  # Next the data frame for the C++-based matrix populator functions
  allstages.list <- theoldpizzle(stageframe, overwrite, repmatrix, finalage = 0, style = 1, cont = 0)
  allstages <- do.call("cbind.data.frame", allstages.list)
  
  maxsize <- max(c(allstages$a.size3, allstages$a.size2n, allstages$a.size2o), na.rm = TRUE)
  
  allstages <- allstages[(which(allstages$b.index321 != -1)),]
  
  # Now we will work up the vital rate models
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
  }
  
  set.seed(randomseed)
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      rvarssummed <- sum(size_proxy$variances[,"vcov"])
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else {
    sizedist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  # Next we create a list of pops, patches, and years in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame(NA, X, as.matrix(year));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- NA
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    listofyears <- cbind.data.frame(NA, NA, as.matrix(year))
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- NA
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  # A few extra tidbits required for the core matrix estimator to work
  total.matrix.dim <- length(stageframe$bin_size_ctr) - 1

  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  # The next line calls the core matrix estimator function
  madsexmadrigal <- lapply(yearlist, jerzeibalowski, allstages, surv_proxy, obs_proxy, size_proxy, repst_proxy,
                           fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jrepst_proxy, inda, indb, indc, 
                           surv_dev, obs_dev, size_dev, repst_dev, fec_dev, jsurv_dev, jobs_dev, jsize_dev, 
                           jrepst_dev, total.matrix.dim, repmod, rvarssummed, sigma, jrvarssummed, jsigma, 
                           maxsize, sizedist, fecdist, negfec)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  if (reduce == TRUE) {
    drops <- .reducer2(a_list, u_list, f_list, ahstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    ahstages <- drops$ahstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA, ahstages = ahstages, 
                 labels = listofyears[,c(1:3)], matrixqc = qcoutput1, modelqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Reduce Matrix Dimensions By Eliminating Empty Stages
#' 
#' \code{.reducer2()} identifies empty stages in a set of ahistorical matrices and
#' removes themfrom all matrices. It is used within \code{\link{flefko2}()} and 
#' \code{\link{rlefko2}()}.
#' 
#' @param A List of population projection matrices, from a \code{lefkoMat} object.
#' @param U List of surviva-transition matrices corresponding to \code{A}.
#' @param F List of fecundity matrices corresponding to \code{A}.
#' @param ahstages Data frame giving the names and identities of ahistorical 
#' stages used to create matrices.
#' 
#' @return Returns a list of reduced \code{A}, \code{U}, and \code{F} matrices, plus the reduced
#' \code{ahstages} object.
#' 
#' @keywords internal
#' @noRd
.reducer2 <- function(A, U, F, ahstages) {
  stagepatterns <- lapply(A, function(X) {
    matrix.sums <- colSums(X) + rowSums(X)
    return(matrix.sums)
  })
  
  used.stages.mat <- do.call("rbind", stagepatterns)
  used.stages.ovr <- colSums(used.stages.mat)
  keep.stages <- which(used.stages.ovr > 0)
  
  Ared <- lapply(A, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Ured <- lapply(U, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Fred <- lapply(F, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  ahstred <- ahstages[keep.stages,]
  
  return(list(A = Ared, U = Ured, F = Fred, ahstages = ahstred))
}

#' Create Raw Historical Population Projection Matrices
#' 
#' \code{rlefko3()} returns raw historical MPMs, including the associated component
#' transition and fecundity matrices, data frames describing the ahistorical stages
#' used and the historical paired stages, and a data frame describing the
#' population, patch, and year associated with each matrix.
#' 
#' @param data A vertical demographic data frame, with variables corresponding to
#' the naming conventions in \code{\link{verticalize3}()}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status of
#' each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param year A variable corresponding to year or observation time, or a set of
#' such values, given in values associated with the year term used in linear model
#' development. Can also equal \code{all}, in which case matrices will be estimated
#' for all years. Defaults to \code{all}.
#' @param pop A variable designating which populations will have matrices
#' estimated. Should be set to specific population names, or to \code{all} if all 
#' populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{all} if
#' matrices should be estimated for all patches. Defaults to \code{all}.
#' @param censor If TRUE, then data will be removed according to the variable set
#' in \code{censorcol}, such that only data with censor values equal to 1 will
#' remain. Defaults to FALSE.
#' @param stages An optional vector denoting the names of the variables within the
#' main vertical dataset coding for the stages of each individual in times \emph{t}+1,
#' \emph{t}, and \emph{t}-1. The names of stages in these variables should match those
#' used in the \code{stageframe} exactly. If left blank, then \code{rlefko3()} will attempt
#' to infer stages by matching values of \code{alive}, \code{size}, \code{repst}, and \code{matst}
#' to characteristics noted in the associated \code{stageframe}.
#' @param alive A vector of names of binomial variables corresponding to status as
#' alive (1) or dead (0) in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' @param size A vector of names of variables coding size in times \emph{t}+1, \emph{t},
#' and \emph{t}-1, respectively. Defaults to \code{c("sizea3", "sizea2", "sizea1")}.
#' @param repst A vector of names of variables coding reproductive status in times
#' \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("repstatus3", 
#' "repstatus2", "repstatus1")}.
#' @param matst A vector of names of variables coding maturity status in times
#' \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("matstatus3", 
#' "matstatus2", "matstatus1")}. Must be supplied if \code{stages} is not provided.
#' @param fec A vector of names of variables coding fecundity in times \emph{t}+1,
#' \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("feca3", "feca2", "feca1")}.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for each
#' potentially new individual (row) born to each reproductive stage (column).
#' Entries act as multipliers on fecundity, with 1 equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function
#' describing transitions to be overwritten either with given values or 
#' with other estimated transitions.
#' @param yearcol The variable name or column number corresponding to time in time
#' \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the identity
#' of the population.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset.
#' @param indivcol The variable name or column number coding individual identity.
#' @param censorcol The variable name or column number denoting the censor status.
#' Only needed if \code{censor = TRUE}.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated exclusively with zero transitions. These are removed only if all row
#' and column sums in ALL matrices estimated equal 0. Defaults to FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return
#' either an object of class \code{lefkoMat}. Output includes:
#'
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and years. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'                           patchidcol = "patch", individcol = "plantid",
#'                           blocksize = 4, sizeacol = "Inf2.04", sizebcol = "Inf.04",
#'                           sizeccol = "Veg.04", repstracol = "Inf.04",
#'                           repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'                           stageassign = cypframe_raw, stagesize = "sizeadded",
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover3r <- overwrite(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL", 
#'                        "SL", "SL", "D", "XSm", "Sm", "D", "XSm", "Sm"), 
#'                        stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", 
#'                        "SL", "SL", "SL", "SL", "SL", "SL", "SL", "SL"),
#'                        stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", 
#'                        "P3", "SL", "P3", "P3", "P3", "SL", "SL", "SL"),
#'                        eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "D", 
#'                        "XSm", "Sm", "D", "XSm", "Sm"), eststage2 = c(NA, NA, NA, 
#'                        NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
#'                        "XSm"), eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'                        "XSm", "XSm", "XSm", "XSm", "XSm", "XSm"), 
#'                        givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, 0.4,
#'                        0.4, NA, NA, NA, NA, NA, NA), type = c("S", "S", "S", 
#'                        "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", 
#'                        "S"))
#' 
#' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, year = "all",
#'                        patch = "all", stages = c("stage3", "stage2", "stage1"),
#'                        size = c("size3added", "size2added", "size1added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover3r,
#'                        yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' summary(cypmatrix3r)
#' 
#' @export
rlefko3 <- function(data, stageframe, year = "all", pop = NA, patch = NA, censor = FALSE, 
                    stages = NA, alive = c("alive3", "alive2", "alive1"), 
                    size = c("sizea3", "sizea2", "sizea1"), repst = c("repstatus3", "repstatus2", "repstatus1"), 
                    matst = c("matstatus3", "matstatus2", "matstatus1"), fec = c("feca3", "feca2", "feca1"), 
                    repmatrix = NA, overwrite = NA, yearcol = NA, popcol = NA, patchcol = NA, indivcol = NA, 
                    censorcol = NA, reduce = FALSE) {
  
  tocensor <- indataset <- alive2 <- popused <- patchused <- yearused <- NULL
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.")
  }
  
  if (!any(class(data) == "hfvdata")) {
    stop("Need original vertical dataset to proceed. This dataset must be in historical vertical format.")
  }
  
  if (all(is.na(stages))) {
    if ((length(alive) != 3)) {
      stop("This function requires stage information for each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
    }
    if ((length(size) != 3)) {
      stop("This function requires stage information for each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
    }
    if (!all(is.na(repst))) {
      if ((length(repst) != 3)) {
        stop("This function requires stage information for each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
      }
    }   
    if (!all(is.na(matst))) {
      if ((length(matst) != 3)) {
        stop("This function requires stage information for each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
      }
    }   
  } else if (length(stages) != 3) {
    stop("This function requires stage information for each of times t+1, t, and t-1.")
  }
  
  if ((length(fec) != 3)) {
    stop("This function requires three variables for fecundity, for each of times t+1, t, and t-1.")
  }
  
  if (any(is.character(year)) & any(class(data) == "hfvdata")) {
    if (is.element(tolower(year), "all")) {
      if (is.character(yearcol)) {
        choicevar <- which(names(data) == yearcol)
        year <- sort(unique(data[,choicevar]))[-1] #Here we remove the 1st year because it has no time t-1 in the dataset
        yearcol <- choicevar
      } else if (all(is.numeric(year))) {
        year <- sort(unique(data[,yearcol]))[-1]
      } else {stop("Cannot understand year designation.")}
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE)) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.")
  }
  
  if (!all(is.element(year, sort(unique(data[,yearcol]))))) {
    stop("Dataset does not contain one or more of the requested years.")
  }
  
  if (all(is.na(repmatrix))) {    #This bit needs to be redone to check for fecundity values in the dataset and use them to determine the proper reproductive stages
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {
      stop("Cannot censor the data without a proper censor variable.")
    }
    
    for (i in c(1:length(censorcol))) {
      if (is.character(censorcol)) {
        data$tocensor <- data[,which(names(data) == censorcol[i])]
      } else {
        data$tocensor <- data[,censorcol[i]]
      }
      data <- subset(data, tocensor == 1)
    }   
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {stop("Need population and patch designation variables to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, sort(unique(patchfolly[,patchcol])));
        names(output) <- c("pop", "patch");
        return(output);
      })
      
      if (length(listofpatches) > 1) {
        listofpatches <- do.call(rbind.data.frame, listofpatches)
      }
    } else {listofpatches <- expand.grid(pop = pops, patch = patch)}
    
    listofyears <- apply(as.matrix(listofpatches), 1, function(X) {
      checkyrdata <- subset(data, popcol = X[1]);
      checkyrdata <- subset(checkyrdata, patchcol = X[2])
      output <- cbind.data.frame(X[1], X[2], sort(unique(checkyrdata[,yearcol]))[-1]);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(patchcol)) {stop("Need patch designation variable to proceed.")}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- sort(unique(data[,patchcol]))
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame(NA, X, sort(unique(checkyrdata[,yearcol]))[-1]);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {stop("Need population designation variable to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, NA, sort(unique(checkyrdata[,yearcol]))[-1]);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    listofyears <- cbind.data.frame(NA, NA, year)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$stage_id <- as.numeric(stageframe$stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  
  data$alive1 <- data[,which(names(data) == alive[3])]
  data$alive2 <- data[,which(names(data) == alive[2])]
  data$alive3 <- data[,which(names(data) == alive[1])]
  
  if (all(is.na(stages))) {
    if (length(size) > 1) {
      data$usedsize1 <- data[,which(names(data) == size[3])]
      data$usedsize2 <- data[,which(names(data) == size[2])]
      data$usedsize3 <- data[,which(names(data) == size[1])]
    } else {
      warning("Without stage columns, lefko3 MPM estimation functions generally require size variables. Failure to include size variables may lead to odd results.")
    }
    if (length(repst) > 1) {
      data$usedrepst1 <- data[,which(names(data) == repst[3])]
      data$usedrepst2 <- data[,which(names(data) == repst[2])]
      data$usedrepst3 <- data[,which(names(data) == repst[1])]
    } 
    if (length(matst) > 1) {
      data$usedmatstatus1 <- data[,which(names(data) == matst[3])]
      data$usedmatstatus2 <- data[,which(names(data) == matst[2])]
      data$usedmatstatus3 <- data[,which(names(data) == matst[1])]
    } 
    
    data$usedstage1 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize1[X])) {
        data$usedsize1[X] <- 0
      }
      mainstages <- intersect(which(stageframe$bin_size_min < data$usedsize1[X]), 
                              which(stageframe$bin_size_max >= data$usedsize1[X]))
      jmstages <- which(stageframe$immstatus == (1 - data$usedmatstatus1[X]))
      obsstages <- which(stageframe$obsstatus == data$obsstatus1[X])
      repstages <- which(stageframe$repstatus == data$repstatus1[X])
      alivestage1 <- which(stageframe$alive == data$alive1[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages)), alivestage1)
      
      if (length(choicestage) == 0) choicestage <- which(stageframe$stage_id == max(stageframe$stage_id))
      
      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      return(as.character(stageframe$stage[choicestage]))
    })
    
    data$usedstage2 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize2[X])) {
        data$usedsize2[X] <- 0
      }
      mainstages <- intersect(which(stageframe$bin_size_min < data$usedsize2[X]), 
                              which(stageframe$bin_size_max >= data$usedsize2[X]))
      jmstages <- which(stageframe$immstatus == (1 - data$usedmatstatus2[X]))
      obsstages <- which(stageframe$obsstatus == data$obsstatus2[X])
      repstages <- which(stageframe$repstatus == data$repstatus2[X])
      
      choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
      
      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      return(as.character(stageframe$stage[choicestage]))
    })
    
    data$usedstage3 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize3[X])) {
        data$usedsize3[X] <- 0
      }
      mainstages <- intersect(which(stageframe$bin_size_min < data$usedsize3[X]), 
                              which(stageframe$bin_size_max >= data$usedsize3[X]))
      jmstages <- which(stageframe$immstatus == (1 - data$usedmatstatus3[X]))
      obsstages <- which(stageframe$obsstatus == data$obsstatus3[X])
      repstages <- which(stageframe$repstatus == data$repstatus3[X])
      alivestage3 <- which(stageframe$alive == data$alive3[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages)), alivestage3)
      
      if (length(choicestage) == 0) choicestage <- which(stageframe$stage_id == max(stageframe$stage_id))
      
      return(as.character(stageframe$stage[choicestage]))
    })
    
  } else if (length(stages) > 1) {
    
    if (is.numeric(stages[2])) {
      data$usedstage1 <- data[, stages[3]]
      data$usedstage2 <- data[, stages[2]]
      data$usedstage3 <- data[, stages[1]]
      
      data$usedstage1 <- as.character(data$usedstage1)
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage1))) {
        data$usedstage1[which(data$usedstage1 == "NotAlive")] <- "Dead"
      }
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    } else {
      data$usedstage1 <- data[,which(names(data) == stages[3])]
      data$usedstage2 <- data[,which(names(data) == stages[2])]
      data$usedstage3 <- data[,which(names(data) == stages[1])]
      
      data$usedstage1 <- as.character(data$usedstage1)
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage1))) {
        data$usedstage1[which(data$usedstage1 == "NotAlive")] <- "Dead"
      }
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    }
    stages.used <- sort(unique(c(data$usedstage2, data$usedstage3)))
    
    if (length(setdiff(stages.used, stageframe$stage)) > 0 & !is.element("NotAlive", stages.used)) {
      stop("Some stages in dataset do not match those detailed in the input stageframe.", .call = FALSE)
    }
  }
  
  if (length(fec) > 1) {
    data$usedfec1 <- data[,which(names(data) == fec[3])]
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    
    data$usedfec1[which(is.na(data$usedfec1))] <- 0
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Lefko3 MPM estimation functions generally require fecundity variables. Failure to include fecundity variables leads to matrices composed only of survival transitions.")
  } 
  
  overwrite <- .overwrite_reassess(overwrite, stageframe, historical = TRUE)
  
  # This section creates stageexpansion9, which is a data frame that holds values for stage transitions from paired stages
  # in times t and t-1 to paired stages in times t and t+1
  majortrial <- theoldpizzle(stageframe, overwrite, repmatrix, 0, 0, 0)
  stageexpansion9 <- do.call("cbind.data.frame", majortrial)
  
  # Stageexpansion3 is a dataframe created to hold values for paired stages in times t and t-1 only
  stageexpansion3 <- cbind.data.frame(expand.grid(size3 = stageframe$bin_size_ctr, size2n = stageframe$bin_size_ctr), 
                                      expand.grid(rep3 = stageframe$repstatus, rep2n = stageframe$repstatus),
                                      expand.grid(indata3 = stageframe$indataset, indata2n = stageframe$indataset),
                                      expand.grid(stage3 = stageframe$stageno, stage2n = stageframe$stageno),
                                      fec32n = c(cbind(rbind(repmatrix, 0), 0)))
  
  stageexpansion3$indata32n <- stageexpansion3$indata3 * stageexpansion3$indata2n
  
  instages <- length(stageframe$stage_id)
  
  stageexpansion3$index21 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    ((stageexpansion3$stage3[X] - 1) + ((stageexpansion3$stage2n[X] - 1) * instages))
  })
  
  stageexpansion3$stcod3 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    stageframe$stage[which(stageframe$stageno == stageexpansion3$stage3[X])]
  })
  stageexpansion3$stcod2 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    stageframe$stage[which(stageframe$stageno == stageexpansion3$stage2n[X])]
  })
  
  # Now we will add a number of indices to the dataset
  data <- subset(data, alive2 == 1)
  
  data$index1 <- apply(as.matrix(data$usedstage1), 1, function(X) {
    stageframe$stageno[which(stageframe$stage == X)]
    
  })
  data$index2 <- apply(as.matrix(data$usedstage2), 1, function(X) {
    stageframe$stageno[which(stageframe$stage == X)]
  })
  data$index3 <- apply(as.matrix(data$usedstage3), 1, function(X) {
    stageframe$stageno[which(stageframe$stage == X)]
  })
  
  data$index321 <- apply(as.matrix(c(1:length(data$usedstage1))), 1, function(X) {
    (data$index3[X] + (data$index2[X] * instages) + (data$index1[X] * instages * instages))
  })
  data$pairindex21 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
    ((data$index2[X] - 1) + ((data$index1[X] - 1) * instages))
  })
  
  data$usedfec2[which(is.na(data$usedfec2))] <- 0
  
  if(is.element(0, unique(data$index1))) {
    warning("Data (stage at time t-1) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  if(is.element(0, unique(data$index2))) {
    warning("Data (stage at time t) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  if(is.element(0, unique(data$index3))) {
    warning("Data (stage at time t+1) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  
  madsexmadrigal <- lapply(yearlist, function(X) {
    passed_data <- data
    if (!is.na(X$pop[1])) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X$pop[1]);
    }
    if (!is.na(X$patch[1])) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X$patch[1]);
    }
    if (!is.na(X$year2[1])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X$year2[1]);
    }
    
    specialpatrolgroup(sge9l = stageexpansion9, sge3 = stageexpansion3, MainData = passed_data, StageFrame = stageframe)
  })
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  indivs <- NA
  if (!all(is.na(indivcol))) {
    if (all(is.character(indivcol))) {indivcol <- which(names(data) == indivcol)[1]}
    indivs <- length(unique(data[,indivcol]))
  }
  qcoutput2 <- c(indivs, dim(data)[1])
  
  morebitstolose <- unique(c(which(stageexpansion3$stage3 == dim(stageframe)[1]), which(stageexpansion3$stage2n == dim(stageframe)[1])))
  stageexpansion3 <- stageexpansion3[-morebitstolose,]
  
  hstages <- stageexpansion3[,c("stage3", "stage2n", "stcod3", "stcod2")]
  names(hstages) = c("stage_id_2", "stage_id_1", "stage_2", "stage_1")
  
  if (reduce == TRUE) {
    drops <- .reducer3(a_list, u_list, f_list, hstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    hstages <- drops$hstages
  }
  
  rownames(hstages) <- c(1:dim(hstages)[1])
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages, 
                 ahstages = stageframe[1:(dim(stageframe)[1] - 1),], labels = listofyears,
                 matrixqc = qcoutput1, dataqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create Raw Ahistorical Population Projection Matrices
#'
#' \code{rlefko2()} returns raw ahistorical MPMs, including the associated
#' componenttransition and fecundity matrices, a data frame describing the
#' ahistorical stages used, and a data frame describing the population, patch, and
#' year associated with each matrix.
#' 
#' @param data A vertical demographic data frame, with variables corresponding 
#' to the naming conventions in \code{\link{verticalize3}()}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param pop A variable designating which populations will have matrices
#' estimated. Should be set to specific population names, or to \code{all} if all 
#' populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{all} if
#' matrices should be estimated for all patches. Defaults to \code{all}.
#' @param censor If TRUE, then data will be removed according to the variable set
#' in \code{censorcol}, such that only data with censor values equal to 1 will
#' remain. Defaults to FALSE.
#' @param stages An optional vector denoting the names of the variables within the
#' main vertical dataset coding for the stages of each individual in times \emph{t}+1,
#' \emph{t}, and \emph{t}-1. The names of stages in these variables should match those
#' used in the \code{stageframe} exactly. If left blank, then \code{rlefko3()} will attempt
#' to infer stages by matching values of \code{alive}, \code{size}, \code{repst}, and \code{matst}
#' to characteristics noted in the associated \code{stageframe}.
#' @param alive A vector of names of binomial variables corresponding to status as
#' alive (1) or dead (0) in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' @param size A vector of names of variables coding size in times \emph{t}+1 and \emph{t},
#' respectively. Defaults to \code{c("sizea3", "sizea2")}.
#' @param repst A vector of names of variables coding reproductive status in times
#' \emph{t}+1 and \emph{t}, respectively. Defaults to \code{c("repstatus3", 
#' "repstatus2")}.
#' @param matst A vector of names of variables coding maturity status in times
#' \emph{t}+1 and \emph{t}, respectively. Defaults to \code{c("matstatus3", 
#' "matstatus2")}. Must be supplied if \code{stages} is not provided.
#' @param fec A vector of names of variables coding fecundity in times \emph{t}+1
#' and \emph{t}, respectively. Defaults to \code{c("feca3", "feca2")}.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for each
#' potentially new individual (row) born to each reproductive stage (column).
#' Entries act as multipliers on fecundity, with 1 equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function
#' describing transitions to be overwritten either with given values or 
#' with other estimated transitions.
#' @param yearcol The variable name or column number corresponding to time in time
#' \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the identity
#' of the population.
#' @param patchcol The variable name or column number corresponding to patch in the
#' dataset.
#' @param indivcol The variable name or column number coding individual identity.
#' @param censorcol The variable name or column number denoting the censor status.
#' Only needed if \code{censor = TRUE}.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated with only zero transitions. These are removed only if all row and
#' column sums in ALL matrices estimated equal 0. Defaults to FALSE.
#' 
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}. This includes:
#' 
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and years. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{Null for ahistorical matrices.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'                           patchidcol = "patch", individcol = "plantid",
#'                           blocksize = 4, sizeacol = "Inf2.04", sizebcol = "Inf.04",
#'                           sizeccol = "Veg.04", repstracol = "Inf.04",
#'                           repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'                           stageassign = cypframe_raw, stagesize = "sizeadded",
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'                        "XSm", "Sm"), stage2 = c("SD", "SD", "P1", "P2", "P3", 
#'                        "SL", "SL", "SL", "SL"), eststage3 = c(NA, NA, NA, NA, 
#'                        NA, NA, "D", "XSm", "Sm"), eststage2 = c(NA, NA, NA, NA, 
#'                        NA, NA, "XSm", "XSm", "XSm"), givenrate = c(0.1, 0.2, 
#'                        0.2, 0.2, 0.25, 0.4, NA, NA, NA), type = c("S", "S", "S",
#'                        "S", "S", "S", "S", "S", "S"))
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, year = "all",
#'                        patch = "all", stages = c("stage3", "stage2", "stage1"),
#'                        size = c("size3added", "size2added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover2r,
#'                        yearcol = "year2", patchcol = "patchid",
#'                        indivcol = "individ")
#'                        
#' cypmatrix2r$A[[1]]
#' 
#' @export
rlefko2 <- function(data, stageframe, year = "all", pop = NA, patch = NA, censor = FALSE, 
                    stages = NA, alive = c("alive3", "alive2"), size = c("sizea3", "sizea2"), 
                    repst = c("repstatus3", "repstatus2"), matst = c("matstatus3", "matstatus2"),
                    fec = c("feca3", "feca2"), repmatrix = NA, overwrite = NA, yearcol = NA, popcol = NA, 
                    patchcol = NA, indivcol = NA, censorcol = NA, reduce = FALSE) {
  
  tocensor <- indataset <- alive2 <- popused <- patchused <- yearused <- NULL
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.")
  }
  
  if (!any(class(data) == "hfvdata")) {
    stop("Need original vertical dataset to proceed. This dataset must be in historical vertical format.")
  }
  
  if (all(is.na(stages))) {
    if (!(length(alive) > 1)) {
      stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
    }
    if (!(length(size) > 1)) {
      stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
    }
    if (!all(is.na(repst))) {
      if (!(length(repst) > 1)) {
        stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
      }
    }   
    if (!all(is.na(matst))) {
      if (!(length(matst) > 1)) {
        stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
      }
    }   
  }
  
  if (!(length(fec) > 1)) {
    stop("This function requires two variables for fecundity, for each of times t+1 and t.")
  }
  
  if (any(is.character(year)) & any(class(data) == "hfvdata")) {
    if (is.element(tolower(year), "all")) {
      if (is.character(yearcol)) {choicevar <- which(names(data) == yearcol);
      year <- sort(unique(data[,choicevar]))
      yearcol <- choicevar
      } else if (all(is.numeric(year))) {
        year <- sort(unique(data[,yearcol]))
      } else {stop("Cannot understand year designation.")}
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE)) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.")
  }
  
  if (!all(is.element(year, sort(unique(data[,yearcol]))))) {
    stop("Dataset does not contain one or more of the requested years.")
  }
  
  if (all(is.na(repmatrix))) {    
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {stop("Cannot censor the data without a proper censor variable.")}
    
    for (i in c(1:length(censorcol))) {
      if (is.character(censorcol)) {
        data$tocensor <- data[,which(names(data) == censorcol[i])]} else {data$tocensor <- data[,censorcol[i]]
        }
      data <- subset(data, tocensor == 1)
    }   
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {stop("Need population and patch designation variables to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, unique(patchfolly[,yearcol]));
        names(output) <- c("pop", "patch");
        return(output);
      })
      
      if (length(listofpatches) > 1) {
        listofpatches <- do.call(rbind.data.frame, listofpatches)
      }
    } else {listofpatches <- expand.grid(pop = pops, patch = patch)}
    
    listofyears <- apply(as.matrix(listofpatches), 1, function(X) {
      checkyrdata <- subset(data, popcol = X[1]);
      checkyrdata <- subset(checkyrdata, patchcol = X[2])
      output <- cbind.data.frame(X[1], X[2], unique(checkyrdata[,yearcol]));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(patchcol)) {stop("Need patch designation variable to proceed.")}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- unique(data[,patchcol])
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame(NA, X, unique(checkyrdata[,yearcol]));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {stop("Need population designation variable to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, NA, unique(checkyrdata[,yearcol]));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    listofyears <- cbind.data.frame(NA, NA, year)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$stage_id <- as.numeric(stageframe$stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  
  data$alive2 <- data[,which(names(data) == alive[2])]
  data$alive3 <- data[,which(names(data) == alive[1])]
  
  instageframe <- subset(stageframe, indataset == 1)
  instages <- dim(stageframe)[1]
  
  if (all(is.na(stages))) {
    if (length(size) > 1) {
      data$usedsize2 <- data[,which(names(data) == size[2])]
      data$usedsize3 <- data[,which(names(data) == size[1])]
    } else {
      warning("Without stage columns, lefko3 MPM estimation functions generally require size variables. Failure to include size variables may lead to odd results.")
    }
    if (length(repst) > 1) {
      data$usedrepst2 <- data[,which(names(data) == repst[2])]
      data$usedrepst3 <- data[,which(names(data) == repst[1])]
    } 
    if (length(matst) > 1) {
      data$usedmatstatus2 <- data[,which(names(data) == matst[2])]
      data$usedmatstatus3 <- data[,which(names(data) == matst[1])]
    } 
    
    data$usedstage2 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize2[X])) {
        data$usedsize2[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize2[X]), 
                              which(instageframe$bin_size_max >= data$usedsize2[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus2[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus2[X])
      repstages <- which(instageframe$repstatus == data$repstatus2[X])
      
      choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
      
      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      return(as.character(instageframe$stage[choicestage]))
    })
    
    data$usedstage3 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize3[X])) {
        data$usedsize3[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize3[X]), 
                              which(instageframe$bin_size_max >= data$usedsize3[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus3[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus3[X])
      repstages <- which(instageframe$repstatus == data$repstatus3[X])
      alivestage3 <- which(instageframe$alive == data$alive3[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages)), alivestage3)
      
      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      if (length(choicestage) == 0) choicestage <- which(instageframe$stage_id == max(instageframe$stage_id))
      
      return(as.character(instageframe$stage[choicestage]))
    })
    
  } else if (length(stages) > 1) {
    if (is.numeric(stages[2])) {
      data$usedstage2 <- data[, stages[2]]
      data$usedstage3 <- data[, stages[1]]
      
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    } else {
      data$usedstage2 <- data[,which(names(data) == stages[2])]
      data$usedstage3 <- data[,which(names(data) == stages[1])]
      
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    }
    stages.used <- sort(unique(c(data$usedstage2, data$usedstage3)))
    
    if (length(setdiff(stages.used, stageframe$stage)) > 0) {
      stop("Some stages in dataset do not match those detailed in the input stageframe.", .call = FALSE)
    }
  }
  
  if (length(fec) > 1) {
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Lefko3 MPM estimation functions generally require fecundity variables. Failure to include fecundity variables leads to matrices composed only of survival transitions.")
  } 
  
  overwrite <- .overwrite_reassess(overwrite, stageframe, historical = FALSE)
  
  # This section creates stageexpansion3, which is a data frame that holds values for stage transitions from time t to t+1
  majortrial <- theoldpizzle(stageframe, overwrite, repmatrix, 0, 1, 0)
  stageexpansion3 <- do.call("cbind.data.frame", majortrial)
  
  # Stageexpansion2 is a dataframe created to hold values for stages in time t only
  stageexpansion2 <- cbind.data.frame(stage2 = as.numeric(stageframe$stageno), size2 = as.numeric(stageframe$bin_size_ctr), 
                                      rep2 = as.numeric(stageframe$repstatus), indata2 = as.numeric(stageframe$indataset),
                                      index2 = (as.numeric(stageframe$stageno) - 1), fec3 = c(rowSums(repmatrix), 0))
  stageexpansion2$fec3[which(stageexpansion2$fec3 > 0)] <- 1
  
  data <- subset(data, alive2 == 1)
  
  data$index2 <- apply(as.matrix(data$usedstage2), 1, function(X) {
    instageframe$stageno[which(instageframe$stage == X)] - 1
  })
  data$index2[which(is.na(data$index2))] <- 0
  data$index3 <- apply(as.matrix(data$usedstage3), 1, function(X) {
    instageframe$stageno[which(instageframe$stage == X)] - 1
  })
  data$index3[which(is.na(data$index3))] <- 0
  data$index32 <- apply(as.matrix(c(1:length(data$usedstage2))), 1, function(X) {
    (data$index3[X] + (data$index2[X] * instages))
  })
  
  if(is.element(0, unique(data$index2))) {
    warning("Data (stage at time t) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  if(is.element(0, unique(data$index3))) {
    warning("Data (stage at time t+1) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  
  # This section runs the core matrix estimator
  madsexmadrigal <- lapply(yearlist, function(X) {
    passed_data <- data
    if (!is.na(X$pop[1])) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X$pop[1]);
    }
    if (!is.na(X$patch[1])) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X$patch[1]);
    }
    if (!is.na(X$year2[1])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X$year2[1]);
    }
    normalpatrolgroup(sge3 = stageexpansion3, sge2 = stageexpansion2, MainData = passed_data, StageFrame = stageframe)
  })
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  indivs <- NA
  if (!all(is.na(indivcol))) {
    if (all(is.character(indivcol))) {indivcol <- which(names(data) == indivcol)[1]}
    indivs <- length(unique(data[,indivcol]))
  }
  qcoutput2 <- c(indivs, dim(data)[1])
  
  ahstages <-  stageframe[1:(dim(stageframe)[1] - 1),]
  
  if (reduce == TRUE) {
    drops <- .reducer2(a_list, u_list, f_list, ahstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    ahstages <- drops$ahstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA, ahstages = ahstages, 
                 labels = listofyears, matrixqc = qcoutput1, dataqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create Function-based Ahistorical Age x Stage Matrix Projection Model
#'
#' Function \code{aflefko2()} returns ahistorical age x stage MPMs corresponding to
#' the patches and years given, including the associated component transition and
#' fecundity matrices, a data frame detailing the characteristics of ahistorical
#' stages, and a data frame characterizing the patch and year combinations
#' corresponding to these matrices. Unlike \code{\link{rlefko2}()} and \code{\link{rlefko3}()},
#' this function currently does not currently distinguish populations.
#'
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{all} if
#' matrices should be estimated for all patches. Defaults to \code{all}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for each
#' potentially new individual (row) born to each reproductive stage (column).
#' Entries act as multipliers on fecundity, with 1 equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function
#' describing transitions to be overwritten either with given values or 
#' with other estimated transitions.
#' @param data The original historical demographic data frame used to estimate
#' vital rates (class \code{hfvdata}). The original data frame is required in order
#' to initialize years and patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate models.
#' If given, then \code{surv_model}, \code{obs_model}, \code{size_model}, \code{repst_model}, \code{fec_model},
#' \code{jsurv_model}, \code{jobs_model}, \code{jsize_model}, \code{jrepst_model}, \code{paramnames},
#' \code{yearcol}, and \code{patchcol} are not required. No models should include size or
#' reproductive status in time \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can be a
#' model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any survival 
#' probability model given in \code{modelsuite}. This model must have been developed 
#' in a modeling exercise testing only the impacts of time \emph{t}.
#' @param obs_model A linear model predicting sprouting or observation probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any observation probability model given in \code{modelsuite}. This model must
#' have been developed in a modeling exercise testing only the impacts of time \emph{t}.
#' @param size_model A linear model predicting size. This can be a model of class
#' \code{glm} or \code{glmer}, both of which require a predicted poisson variable under a 
#' log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian response is 
#' assumed. If given, then will overwrite any size model given in \code{modelsuite}.  
#' This model must have been developed in a modeling exercise testing only the
#' impacts of time \emph{t}.
#' @param repst_model A linear model predicting reproduction probability. This can
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any reproduction 
#' probability model given in \code{modelsuite}. This model must have been developed 
#' in a modeling exercise testing only the impacts of time \emph{t}.
#' @param fec_model A linear model predicting fecundity. This can be a model of
#' class \code{glm} or \code{glmer}, and requires a predicted poisson variable under a 
#' log link. If given, then will overwrite any fecundity model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing only the
#' impacts of time \emph{t}.
#' @param jsurv_model A linear model predicting juvenile survival probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a predicted 
#' binomial variable under a logit link. If given, then will overwrite any 
#' juvenile survival probability model given in \code{modelsuite}. This model must 
#' have been developed in a modeling exercise testing only the impacts of time \emph{t}.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any juvenile observation probability model given in \code{modelsuite}. This model
#' must have been developed in a modeling exercise testing only the impacts of 
#' time \emph{t}.
#' @param jsize_model A linear model predicting juvenile size. This can be a model
#' of class \code{glm} or \code{glmer}, both of which require a predicted poisson variable 
#' under a log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian 
#' response is assumed. If given, then will overwrite any juvenile size model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing only the impacts of time \emph{t}.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in the previous year. This can be a model 
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial variable under a 
#' logit link. If given, then will overwrite any reproduction probability model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing only the impacts of time \emph{t}.
#' @param paramnames A dataframe with two columns, the first showing the general
#' model terms that will be used in matrix creation, and the second showing the
#' equivalent terms used in modeling. Only required if \code{modelsuite} is not 
#' supplied.
#' @param inda A numeric value to use for individual covariate a. Defaults to 0.
#' @param indb A numeric value to use for individual covariate b. Defaults to 0.
#' @param indc A numeric value to use for individual covariate c. Defaults to 0.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile reproduction probability.
#' @param repmod A scalar multiplier of fecundity. Defaults to 1.
#' @param yearcol The variable name or column number corresponding to year in time
#' \emph{t} in the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for missing
#' patches within vital rate models should be estimated as random intercepts.
#' Defaults to FALSE, in which case missing time step coefficients are set to 0.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing patch coefficients are set
#' to 0.
#' @param final_age The final age to model in the matrix, where the first age will
#' be age 0.
#' @param continue A logical value designating whether to allow continued survival
#' of individuals going past the final age, using the demographic characteristics
#' of the final age.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing time step and patch coefficients, if either \code{year.as.random} or
#' \code{patch.as.random} is set to TRUE. Defaults to 0.
#' @param negfec A logical value denoting whether fecundity values estimated to be
#' negative should be reset to 0. Defaults to FALSE.
#' @param reduce A logical value denoting whether to remove ahistorical stages
#' associated solely with 0 transitions. These are only removed in cases where
#' the associated row and column sums in ALL matrices estimated equal 0. 
#' Defaults to FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}. Output includes:
#'
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' years.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}.}
#' \item{hstages}{Null for ahistorical matrices.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the patch and year of each matrix in order.
#' In \code{aflefko2()}, only one population may be analyzed at once, and so \code{pop = NA}}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the modelsuite input.}
#' 
#' Please note that this function will yield incorrect estimates if the models
#' utilized incorporate state in time \emph{t}-1. Only use models developed testing
#' ahistorical effects.
#'
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", "Sz5nr",
#'                  "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", "Sz4r",
#'                  "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' minima <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' maxima <- c(NA, 1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                         obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                         indataset = indataset, binhalfwidth = binvec, propstatus = propvector,
#'                         minage = minima, maxage = maxima)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                            individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                            sizeacol = "lnVol88", repstracol = "Intactseed88",
#'                            fecacol = "Intactseed88", deadacol = "Dead1988",
#'                            nonobsacol = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathrepmln <- matrix(0, 21, 21)
#' lathrepmln[1, c(13:21)] <- 0.345
#' lathrepmln[2, c(13:21)] <- 0.054
#' 
#' lathover2 <- overwrite(stage3 = c("Sd", "Sdl"), stage2 = c("Sd", "Sd"),
#'                        givenrate = c(0.345, 0.054))
#' 
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "lme4", suite = "main", 
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl",bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", age = "obsage", year.as.random = TRUE, 
#'                              patch.as.random = TRUE, show.model.tables = TRUE, quiet = TRUE)
#'                              
#' lathmat2age <- aflefko2(year = "all", patch = "all", stageframe = lathframeln, 
#'                         modelsuite = lathmodelsln2, data = lathvertln, 
#'                         repmatrix = lathrepmln, overwrite = lathover2,
#'                         patchcol = "patchid", yearcol = "year2", year.as.random = FALSE,
#'                         patch.as.random = FALSE, final_age = 2, continue = TRUE, reduce = FALSE)
#' 
#' summary(lathmat2age)
#' }
#' 
#' @export
#' 
aflefko2 <- function(year = 1, patch = NA, stageframe, repmatrix = NA, overwrite = NA, 
                     data = NA, modelsuite = NA, surv_model = NA, obs_model = NA, size_model = NA, 
                     repst_model = NA, fec_model = NA, jsurv_model = NA, jobs_model = NA, 
                     jsize_model = NA, jrepst_model = NA, paramnames = NA, inda = 0, indb = 0, 
                     indc = 0, surv_dev = 0, obs_dev = 0, size_dev = 0, repst_dev = 0, fec_dev = 0, 
                     jsurv_dev = 0, jobs_dev = 0, jsize_dev = 0, jrepst_dev = 0, repmod = 1, 
                     yearcol = "year2", patchcol = "patchid", year.as.random = FALSE, patch.as.random = FALSE, 
                     final_age = 10, continue = TRUE, randomseed = 0, negfec = FALSE, reduce = FALSE) {
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warnings("Function may not work properly without a dataframe of model parameters or equivalents supplied either through modelsuite or through the paramnames input parameter.")
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {stop("Need original vertical dataset to set proper limits on year and patch.", call. = FALSE)}
  if (!any(class(data) == "data.frame")) {stop("Need original vertical dataset used in modeling to proceed.", call. = FALSE)}
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element(tolower(year), "all")) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.", call. = FALSE)
  }
  
  if (is.na(patch) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch option when using a modelsuite in which patch is designated.")
  }
  
  if (is.character(patchcol)) {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element(tolower(patch), "all")) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (all(is.na(repmatrix))) {    #This bit needs to be redone to check for fecundity values in the dataset and use them to determine the proper reproductive stages
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (!is.numeric(inda) | !is.numeric(indb) | !is.numeric(indc)) {
    stop("Individual covariate values must equal valid numbers.", call. = FALSE)
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$size)))))) {
    stop("Function aflefko2() requires size to be numeric rather than categorical.")
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$stage_id <- as.numeric(stageframe$stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$alive <- as.numeric(stageframe$alive)
  stageframe$min_age <- as.numeric(stageframe$min_age)
  stageframe$max_age <- as.numeric(stageframe$max_age)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  stageframe$fullstage <- apply(as.matrix(c(1:dim(stageframe)[1])), 1, function(X) {
    paste(stageframe$bin_size_ctr[X], stageframe$repstatus[X])
  })
  
  instages <- length(stageframe$stage_id)
  
  stageframe$min_age[which(is.na(stageframe$min_age))] <- 0
  
  if (final_age < max(stageframe$max_age, na.rm = TRUE)) {
    warning("Value of final_age will be adjusted to equal the maximum max_age value in the stageframe supplied.")
    final_age <- max(stageframe$max_age, na.rm = TRUE)
  }
  
  # Next the data frame for the C++-based matrix populator functions
  allstages.list <- theoldpizzle(stageframe, overwrite, repmatrix, finalage = final_age, 
                                 style = 2, cont = continue)
  allstages <- do.call("cbind.data.frame", allstages.list)
  
  maxsize <- max(c(allstages$a.size3, allstages$a.size2n), na.rm = TRUE)
  
  allstages <- allstages[(which(allstages$c.aliveandequal != -1)),]
  
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
  }
  
  set.seed(randomseed)
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      rvarssummed <- sum(size_proxy$variances[,"vcov"])
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else {
    sizedist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random)
  
  # This creates a list of pop, patch, and year in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame(NA, X, as.matrix(year));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- NA
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    
    listofyears <- cbind.data.frame(NA, NA, as.matrix(year))
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- NA
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  total.matrix.dim <- (instages - 1) * (final_age + 1)

  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  madsexmadrigal <- lapply(yearlist, jerzeibalowski, allstages, surv_proxy, obs_proxy, size_proxy, repst_proxy,
                           fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jrepst_proxy, inda, indb, indc, 
                           surv_dev, obs_dev, size_dev, repst_dev, fec_dev, jsurv_dev, jobs_dev, jsize_dev, 
                           jrepst_dev, total.matrix.dim, repmod, rvarssummed, sigma, jrvarssummed, jsigma, 
                           maxsize, sizedist, fecdist, negfec)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  if (reduce == TRUE) {
    drops <- .reducer2(a_list, u_list, f_list, ahstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    ahstages <- drops$ahstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA, ahstages = ahstages, 
                 labels = listofyears[,c(1:3)], matrixqc = qcoutput1, modelqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Summary of Class "lefkoMat"
#'
#' A function to simplify the viewing of basic information describing the matrices
#' produced through functions \code{\link{flefko3}()}, \code{\link{flefko2}()}, \code{\link{rlefko3}()},
#' and \code{\link{rlefko2}()}.
#' 
#' @param object An object of class \code{lefkoMat}.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the number of each type of matrix, the
#' number of annual matrices, the number of estimated (non-zero) elements across
#' all matrices and per matrix, the number of unique transitions in the dataset,
#' and the number of individuals.
#' 
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'                           patchidcol = "patch", individcol = "plantid",
#'                           blocksize = 4, sizeacol = "Inf2.04", sizebcol = "Inf.04",
#'                           sizeccol = "Veg.04", repstracol = "Inf.04",
#'                           repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'                           stageassign = cypframe_raw, stagesize = "sizeadded",
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'                        "XSm", "Sm"), stage2 = c("SD", "SD", "P1", "P2", "P3", 
#'                        "SL", "SL", "SL", "SL"), eststage3 = c(NA, NA, NA, NA, 
#'                        NA, NA, "D", "XSm", "Sm"), eststage2 = c(NA, NA, NA, NA, 
#'                        NA, NA, "XSm", "XSm", "XSm"), givenrate = c(0.1, 0.2, 
#'                        0.2, 0.2, 0.25, 0.4, NA, NA, NA), type = c("S", "S", "S",
#'                        "S", "S", "S", "S", "S", "S"))
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, year = "all",
#'                        patch = "all", stages = c("stage3", "stage2", "stage1"),
#'                        size = c("size3added", "size2added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover2r,
#'                        yearcol = "year2", patchcol = "patchid",
#'                        indivcol = "individ")
#'                        
#' summary(cypmatrix2r)
#' 
#' @export
summary.lefkoMat <- function(object, ...) {
  
  matrices <- object
  
  matdim <- dim(matrices$A[[1]])[1]
  
  mqca <- matrices$matrixqc[1]
  mqcb <- matrices$matrixqc[2]
  mqcc <- matrices$matrixqc[3]
  
  if (!all(is.na(matrices$hstages))) {
    histmark <- "historical"
  } else {
    histmark <- "ahistorical"
  }
  
  mqac <- mqca / mqcc
  if (mqac != floor(mqac)) mqac <- round(mqac, digits = 3)

  mqbc <- mqcb / mqcc
  if (mqbc != floor(mqbc)) mqbc <- round(mqbc, digits = 3)
  
  writeLines(paste0("\nThis ", histmark, " lefkoMat object contains ", mqcc, " matrices."))
  writeLines(paste0("\nEach matrix is a square matrix with ", matdim, " rows and columns, and a total of ", matdim*matdim, " elements."))
  writeLines(paste0("A total of ", mqca, " survival transitions were estimated, with ", 
                    mqac, " per matrix."))
  writeLines(paste0("A total of ", mqcb, " fecundity transitions were estimated, with ", 
                    mqbc, " per matrix."))
  
  if (is.element("dataqc", names(matrices))) {
    dqca <- matrices$dataqc[1]
    dqcb <- matrices$dataqc[2]
    
    writeLines(paste0("\nThe dataset contains a total of ", dqca, " unique individuals and ", dqcb, " unique transitions."))
  }
  
  if (is.element("modelqc", names(matrices))) {
    moqc12 <- matrices$modelqc[1,2]
    moqc22 <- matrices$modelqc[2,2]
    moqc32 <- matrices$modelqc[3,2]
    moqc42 <- matrices$modelqc[4,2]
    moqc52 <- matrices$modelqc[5,2]
    moqc62 <- matrices$modelqc[6,2]
    moqc72 <- matrices$modelqc[7,2]
    moqc82 <- matrices$modelqc[8,2]
    moqc92 <- matrices$modelqc[9,2]
    
    moqc13 <- matrices$modelqc[1,3]
    moqc23 <- matrices$modelqc[2,3]
    moqc33 <- matrices$modelqc[3,3]
    moqc43 <- matrices$modelqc[4,3]
    moqc53 <- matrices$modelqc[5,3]
    moqc63 <- matrices$modelqc[6,3]
    moqc73 <- matrices$modelqc[7,3]
    moqc83 <- matrices$modelqc[8,3]
    moqc93 <- matrices$modelqc[9,3]
    
    writeLines("\nVital rate modeling quality control:\n")

    if (moqc12 > 0) {
      writeLines(paste0("Survival estimated with ", moqc12, " individuals and ", moqc13, " individual transitions."))
    } else {
      writeLines("Survival not estimated.")
    }
    
    if (moqc22 > 0) {
      writeLines(paste0("Observation estimated with ", moqc22, " individuals and ", moqc23, " individual transitions."))
    } else {
      writeLines("Observation probability not estimated.")
    }
    
    if (moqc32 > 0) {
      writeLines(paste0("Size estimated with ", moqc32, " individuals and ", moqc33, " individual transitions."))
    } else {
      writeLines("Size transition not estimated.")
    }
    
    if (moqc42 > 0) {
      writeLines(paste0("Reproductive status estimated with ", moqc42, " individuals and ", moqc43, " individual transitions."))
    } else {
      writeLines("Reproduction probability not estimated.")
    }
    
    if (moqc52 > 0) {
      writeLines(paste0("Fecundity estimated with ", moqc52, " individuals and ", moqc53, " individual transitions."))
    } else {
      writeLines("Fecundity not estimated.")
    }
    
    if (moqc62 > 0) {
      writeLines(paste0("Juvenile survival estimated with ", moqc62, " individuals and ", moqc63, " individual transitions."))
    } else {
      writeLines("Juvenile survival not estimated.")
    }
    
    if (moqc72 > 0) {
      writeLines(paste0("Juvenile observation estimated with ", moqc72, " individuals and ", moqc73, " individual transitions."))
    } else {
      writeLines("Juvenile observation probability not estimated.")
    }
    
    if (moqc82 > 0) {
      writeLines(paste0("Juvenile size estimated with ", moqc82, " individuals and ", moqc83, " individual transitions."))
    } else {
      writeLines("Juvenile size transition not estimated.")
    }
    
    if (moqc92 > 0) {
      writeLines(paste0("Juvenile reproduction estimated with ", moqc92, " individuals and ", moqc93, " individual transitions."))
    } else {
      writeLines("Juvenile reproduction probability not estimated.")
    }
  }
  return()
}

