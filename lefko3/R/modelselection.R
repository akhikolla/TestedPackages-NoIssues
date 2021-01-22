#' Develop Best-fit Vital Rate Estimation Models For Matrix Development
#' 
#' \code{modelsearch()} returns both a table of vital rate estimating models and a
#' best-fit model for each major vital rate estimated. The final output can be used
#' as input in other functions within this package.
#' 
#' @param data The vertical dataset to be used for analysis. This dataset should 
#' be of class \code{hfvdata}, but can also be a data frame formatted similarly to
#' the output format provided by functions \code{\link{verticalize3}()} or \code{\link{historicalize3}()},
#' as long as all needed variables are properly designated.
#' @param historical A logical variable denoting whether to assess the effects of
#' state in time \emph{t}-1 in addition to state in time \emph{t}. Defaults to TRUE.
#' @param approach The statistical approach to be taken for model building. The 
#' default is \code{lme4}, which uses the mixed model approach utilized in 
#' package \code{lme4}. Other options include \code{glm}, which uses \code{lm}, \code{glm},
#' and related functions in base R.
#' @param suite This describes the global model for each vital rate estimation and
#' has the following possible values: \code{full}, includes main effects and all
#' two-way interactions of size and reproductive status; \code{main}, includes main
#' effects only of size and reproductive status; \code{size}, includes only size
#' (also interactions between size in historical model); \code{rep}, includes only
#' reproductive status (also interactions between status in historical model); 
#' \code{cons}, all vital rates estimated only as y-intercepts. If \code{approach = "glm"}
#' and \code{year.as.random = FALSE}, then year is also included as a fixed effect,
#' and, in the case of \code{full}, included in two-way interactions. Defaults to \code{size}.
#' @param bestfit A variable indicating the model selection criterion for the
#' choice of best-fit model. The default is \code{AICc&k}, which chooses the 
#' best-fit model as the model with the lowest AICc or, if not the same model, then
#' the model that has the lowest degrees of freedom among models with \eqn{\Delta AICc <= 2.0}.
#' Alternatively, \code{AICc} may be chosen, in which case the best-fit model is simply
#' the model with the lowest AICc value.
#' @param vitalrates A vector describing which vital rates will be estimated via
#' linear modeling, with the following options: \code{surv}, survival probability;
#' \code{obs}, observation probability; \code{size}, overall size; \code{repst}, 
#' probability of reproducing; and \code{fec}, amount of reproduction (overall 
#' fecundity). Defaults to \code{c("surv", "size", "fec")}.
#' @param surv A vector indicating the variable names coding for status as alive or
#' dead in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("alive3", "alive2", "alive1")}.
#' @param obs A vector indicating the variable names coding for observation status
#' in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("obsstatus3", "obsstatus2", "obsstatus1")}.
#' @param size A vector indicating the variable names coding for size in times \emph{t}+1,
#' \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("sizea3", "sizea2", "sizea1")}.
#' @param repst A vector indicating the variable names coding for reproductive
#' status in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("repstatus3", "repstatus2", "repstatus1")}.
#' @param fec A vector indicating the variable names coding for fecundity in times 
#' \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("feca3", "feca2", "feca1")}.
#' @param stage A vector indicating the variables coding for stage in times \emph{t}+1, 
#' \emph{t}, and \emph{t}-1. Defaults to \code{c("stage3", "stage2", "stage1")}.
#' @param indiv A variable indicating the variable name coding individual identity.
#' Defaults to \code{individ}.
#' @param patch A variable indicating the variable name coding for patch, where
#' patches are defined as permanent subgroups within the study population. Defaults
#' to NA.
#' @param year A variable indicating the variable coding for observation time in 
#' time \emph{t}. Defaults to \code{year2}.
#' @param sizedist The probability distribution used to model size. Options include
#' \code{gaussian} for the Normal distribution (default), \code{poisson} for the Poisson
#' distribution, and \code{negbin} for the negative binomial distribution.
#' @param fecdist The probability distribution used to model fecundity. Options
#' include \code{gaussian} for the Normal distribution (default), \code{poisson} for
#' the Poisson distribution, and \code{negbin} for the negative binomial distribution.
#' @param size.zero A logical variable indicating whether size distribution 
#' should be zero-inflated. Only applies to Poisson and negative binomial 
#' distributions. Defaults to FALSE.
#' @param fec.zero A logical variable indicating whether fecundity distribution 
#' should be zero-inflated. Only applies to Poisson and negative binomial 
#' distributions. Defaults to FALSE.
#' @param patch.as.random If set to TRUE and \code{approach = "lme4"}, then \code{patch}
#' is included as a random factor. If set to FALSE and \code{approach = "glm"}, then 
#' \code{patch} is included as a fixed factor. All other combinations of logical 
#' value and \code{approach} lead to \code{patch} not being included in modeling.
#' Defaults to TRUE.
#' @param year.as.random If set to TRUE and \code{approach = "lme4"}, then \code{year} is 
#' included as a random factor. If set to FALSE, then \code{year} is included as a 
#' fixed factor. All other combinations of logical value and \code{approach} lead to 
#' \code{year} not being included in modeling. Defaults to TRUE.
#' @param juvestimate An optional variable denoting the stage name of the juvenile
#' stage in the vertical dataset. If not NA, and \code{stage} is also given (see
#' below), then vital rates listed in \code{vitalrates} other than \code{fec} will also
#' be estimated from the juvenile stage to all adult stages. Defaults to NA, in
#' which case juvenile vital rates are not estimated.
#' @param juvsize A logical variable denoting whether size should be used as a term
#' in models involving transition from the juvenile stage. Defaults to FALSE, and
#' is only used if \code{juvestimate} does not equal NA.
#' @param fectime A variable indicating which year of fecundity to use as the
#' response term in fecundity models. Options include \code{2}, which refers to time \emph{t},
#' and \code{3}, which refers to time \emph{t}+1. Defaults to \code{2}.
#' @param censor A vector denoting the names of censoring variables in the dataset,
#' in order from time \emph{t}+1, followed by time \emph{t}, and lastly followed by time \emph{t}-1.
#' Defaults to NA.
#' @param age Designates the name of the variable corresponding to age in the
#' vertical dataset. Defaults to NA, in which case age is not included in linear
#' models. Should only be used if building age x stage matrices.
#' @param indcova Vector designating the names in times \emph{t}+1, \emph{t}, and
#' \emph{t}-1 of an individual covariate. Defaults to NA.
#' @param indcovb Vector designating the names in times \emph{t}+1, \emph{t}, and
#' \emph{t}-1 of an individual covariate. Defaults to NA.
#' @param indcovc Vector designating the names in times \emph{t}+1, \emph{t}, and
#' \emph{t}-1 of an individual covariate. Defaults to NA.
#' @param show.model.tables If set to TRUE, then includes full modeling tables
#' in the output. Defaults to TRUE.
#' @param global.only If set to TRUE, then only global models will be built and
#' evaluated. Defaults to FALSE.
#' @param quiet If set to TRUE, then model building and selection will proceed
#' without warnings and diagnostic messages being issued. Note that this will
#' not affect warnings and messages generated as models themselves are tested.
#' Defaults to FALSE.
#' 
#' @return This function yields an object of class \code{lefkoMod}, which is a list
#' in which the first 9 elements are the best-fit models for survival, observation
#' status, size, reproductive status, fecundity, juvenile survival, juvenile
#' observation, juvenile size, and juvenile transition to reproduction, 
#' respectively, followed by 9 elements corresponding to the model tables for each
#' of these vital rates, in order, followed by a single character element denoting
#' the criterion used for model selection, and ending on a quality control vector:
#' 
#' \item{survival_model}{Best-fit model of the binomial probability of survival
#' from time \emph{t} to time \emph{t}+1. Defaults to 1.}
#' \item{observation_model}{Best-fit model of the binomial probability of 
#' observation in time \emph{t}+1 given survival to that time. Defaults to 1.}
#' \item{size_model}{Best-fit model of size in time \emph{t}+1 given survival to
#' and observation in that time. Defaults to 1.}
#' \item{repstatus_model}{Best-fit model of the binomial probability of
#' reproduction in time \emph{t}+1, given survival to and observation in that time.
#' Defaults to 1.}
#' \item{fecundity_model}{Best-fit model of fecundity in time \emph{t}+1 given 
#' survival to, and observation and reproduction in that time. Defaults to 1.}
#' \item{juv_survival_model}{Best-fit model of the binomial probability of survival
#' from time \emph{t} to time \emph{t}+1 of an immature individual. Defaults to 1.}
#' \item{juv_observation_model}{Best-fit model of the binomial probability of 
#' observation in time \emph{t}+1 given survival to that time of an immature
#' individual. Defaults to 1.}
#' \item{juv_size_model}{Best-fit model of size in time \emph{t}+1 given survival to
#' and observation in that time of an immature individual. Defaults to 1.}
#' \item{juv_reproduction_model}{Best-fit model of the binomial probability of
#' reproduction in time \emph{t}+1, given survival to and observation in that time
#' of an individual that was immature in time \emph{t}. This model is technically
#' not a model of reproduction probability for individuals that are immature,
#' rather reproduction probability here is given for individuals that are mature in
#' time \emph{t}+1 but immature in time \emph{t}. Defaults to 1.}
#' \item{survival_table}{Full dredge model table of survival probability.}
#' \item{observation_table}{Full dredge model table of observationprobability.}
#' \item{size_table}{Full dredge model table of size.}
#' \item{repstatus_table}{Full dredge model table of reproduction probability.}
#' \item{fecundity_table}{Full dredge model table of fecundity.}
#' \item{juv_survival_table}{Full dredge model table of immature survival 
#' probability.}
#' \item{juv_observation_table}{Full dredge model table of immature observation
#' probability.}
#' \item{juv_size_table}{Full dredge model table of immature size.}
#' \item{juv_reproduction_table}{Full dredge model table of immature reproduction
#' probability.}
#' \item{criterion}{Vharacter variable denoting the criterion used to determine
#' the best-fit model.}
#' \item{qc}{Data frame with three variables: 1) Name of vital rate, 2) number of
#' individuals used to model that vital rate, and 3) number of individual
#' transitions used to model that vital rate.}
#' 
#' The mechanics governing model building are fairly robust to errors and
#' exceptions. The function attempts to build global models, and simplifies models
#' automatically should model building fail. Model building proceeds through the 
#' functions \code{\link[stats]{lm}()} (GLM with Gaussian response), \code{\link[stats]{glm}()} (GLM 
#' with Poisson or binomial response), \code{\link[MASS]{glm.nb}()} (GLM with negative binomial 
#' response), \code{\link[pscl]{zeroinfl}()} (zero-inflated Poisson or negative binomial response), 
#' \code{\link[lme4]{lmer}()} (mixed model with Gaussian response), \code{\link[lme4]{glmer}()} 
#' (mixed model with binomial or Poisson response), \code{\link[glmmTMB]{glmmTMB}()} (mixed
#' model with negative binomial, zero-inflated negative binomial, or zero-inflated 
#' Poisson response). See documentation related to these functions for further 
#' information.
#' 
#' Exhaustive model building and selection proceeds via the \code{\link[MuMIn]{dredge}()}
#' function in package \code{MuMIn}. This function is verbose, so that any errors and 
#' warnings developed during model building, model analysis, and model selection
#' can be found and dealt with. Interpretations of errors during global model
#' analysis may be found in documentation in for the functions and packages
#' mentioned. Package \code{MuMIn} is used for model dredging (see \link[MuMIn]{dredge}),
#' and errors and warnings during dredging can be interpreted using the
#' documentation for that package. Errors occurring during dredging lead to the
#' adoption of the global model as the best-fit, and the user should view all
#' logged errors and warnings to determine the best way to proceed. The 
#' \code{quiet = TRUE} option can be used to silence dredge warnings, but users should
#' note that automated model selection can be viewed as a black box, and so great
#' care should be taken to ensure that the models run make biological sense, and
#' that model quality is prioritized.
#' 
#' Exhaustive model selection through dredging works best with larger datasets and
#' fewer tested parameters. Setting \code{suite = "full"} may initiate a dredge
#' that takes a dramatically long time, particularly if the model is historical,
#' individual covariates are used, or a zero-inflated distribution is assumed. In
#' such cases, the number of models built and tested will run at least in the millions.
#' Small datasets will also increase the error associated with these tests, leading to
#' adoption of simpler models overall. We do not yet offer a parallelization option
#' for function \code{modelsearch()}, but plan to offer one in the future to speed
#' this process up for particularly large global models.
#'
#' Care must be taken to build models that test the impacts of state in time \emph{t}-1
#' for historical models, and that do not test these impacts for ahistorical
#' models. Ahistorical matrix modeling particularly will yield biased transition
#' estimates if historical terms from models are ignored. This can be dealt with
#' at the start of modeling by setting \code{historical = FALSE} for the
#' ahistorical case, and \code{historical = TRUE} for the historical case.
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
#'                            sizeacol = "lnVol88", repstracol = "FCODE88",
#'                            fecacol = "Intactseed88", deadacol = "Intactseed88",
#'                            nonobsacol = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "lme4", suite = "main",
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl", bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE)
#' 
#' lathmodelsln2
#' }
#' 
#' @export
modelsearch <- function(data, historical = TRUE, approach = "mixed", suite = "size", bestfit = "AICc&k",
                        vitalrates = c("surv", "size", "fec"), surv = c("alive3", "alive2", "alive1"),
                        obs = c("obsstatus3", "obsstatus2", "obsstatus1"), size = c("sizea3", "sizea2", "sizea1"),
                        repst = c("repstatus3", "repstatus2", "repstatus1"), fec = c("feca3", "feca2", "feca1"),
                        stage = c("stage3", "stage2", "stage1"), indiv = "individ", patch = NA, year = "year2",
                        sizedist = "gaussian", fecdist = "gaussian", size.zero = FALSE, fec.zero = FALSE,
                        patch.as.random = TRUE, year.as.random = TRUE, juvestimate = NA, juvsize = FALSE, 
                        fectime = 2, censor = NA, age = NA, indcova = NA, indcovb = NA, indcovc = NA, 
                        show.model.tables = TRUE, global.only = FALSE, quiet = FALSE) {
  
  old <- options() #This function requires changes to options(na.action) in order for the lme4::dredge routines to work properly
  on.exit(options(old)) #This will reset options() to user originals when the function exits
  
  censor1 <- censor2 <- censor3 <- NULL
  
  #Input testing, input standardization, and exception handling
  if (all(class(data) != "hfvdata")) {warning("This function was made to work with standardized historically-formatted vertical datasets, as provided by the verticalize() and historicalize() functions. Failure to format the input data properly and designate needed variables appropriately may result in nonsensical output.")}
  
  if (!requireNamespace("MuMIn", quietly = TRUE)) {stop("Package MuMIn needed for this function to work. Please install it.", call. = FALSE)}
  if (!requireNamespace("stringr", quietly = TRUE)) {stop("Package stringr needed for this function to work. Please install it.", call. = FALSE)}
  
  approach <- tolower(approach)
  sizedist <- tolower(sizedist)
  fecdist <- tolower(fecdist)
  
  if (approach == "lme4") {approach <- "mixed"}
  
  if (approach == "mixed" & !requireNamespace("lme4", quietly = TRUE)) {
    if (sizedist == "negbin" & !requireNamespace("glmmTMB", quietly = TRUE)) {stop("Package glmmTMB needed to develop mixed size models with a negative binomial distribution.")}
    if (fecdist == "negbin" & !requireNamespace("glmmTMB", quietly = TRUE)) {stop("Package glmmTMB needed to develop mixed fecundity models with a negative binomial distribution.")}
    stop("Package lme4 needed for this function to work. Please install it.", call. = FALSE)
  }
  distoptions <- c("gaussian", "poisson", "negbin")
  packoptions <- c("mixed", "glm") #The mixed option now handles all mixed models
  
  if (!is.element(approach, packoptions)) {stop("Please enter a valid approach, currently either 'mixed' or 'glm'.", call. = FALSE)}
  if (!is.element(sizedist, distoptions)) {stop("Please enter a valid assumed size distribution, currently either gaussian, poisson, or negbin.", call. = FALSE)}
  if (!is.element(fecdist, distoptions)) {stop("Please enter a valid assumed fecundity distribution, currently either gaussian, poisson, or negbin.", call. = FALSE)}
  
  if (length(censor) > 3) {
    stop("Censor variables should be included either as 1 variable per row in the historical data file (1 variable in the dataset), or as 1 variable per timestep within each historical data file (2 or 3 variables in the dataset). No more than 3 variables are allowed, and if more than one are supplied, then they are assumed to be in order of time t+1, time t, and time t-1, respectively.", call. = FALSE)
  }
  if (length(indiv) > 1) {stop("Only one individual identification variable is allowed.", call. = FALSE)}
  if (length(year) > 1) {stop("Only one time variable is allowed, and it must refer to time t.", call. = FALSE)}
  if (length(patch) > 1) {stop("Only one patch variable is allowed.", call. = FALSE)}
  
  if (is.element("surv", vitalrates)) {
    if (length(surv) > 3 | length(surv) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) survival variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("obs", vitalrates)) {
    if (length(obs) > 3 | length(obs) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) observation variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("size", vitalrates)) {
    if (length(size) > 3 | length(size) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) size variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("repst", vitalrates)) {
    if (length(repst) > 3 | length(repst) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) reproductive status variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("fec", vitalrates)) {
    if (length(fec) > 3 | length(fec) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) fecundity variables from the dataset as input parameters.", call. = FALSE)}
  }
  
  if (fectime != 2 & fectime != 3) {
    stop("The fectime option must equal either 2 or 3, depending on whether the fecundity response term is for time t or time t+1, respectively.", call. = FALSE)
  }
  
  if (!is.na(age)) {
    if (length(which(names(data) == age)) == 0) {
      stop("Variable age must equal either NA or the exact name of the variable denoting age in the dataset.", call. = FALSE)
    } else {
      agecol <- which(names(data) == age)
    }
  } else {age <- "none"}
  
  if (any(!is.na(indcova))) {
    if (length(indcova) > 3) {
      warning("Vector indcova holds the names of an individual covariate across times t+1, t, and t-1. Will use only the first three elements.")
      indcova <- indcova[1:3]
    } else if (length(indcova) == 1) {
      warning("Vector indcova requires the names of an individual covariate across times t+1, t, and t-1. Since only 1 name was supplied, will assume that this name corresponds to status in time t+1 and not use the variable.")
      indcova <- c("none", "none", "none")
    } else {
      if (length(which(names(data) == indcova[2])) == 0 && indcova[2] != "none") {
        stop("Variable indcova must equal either NA or a vector of up to 3 exact names of an individual covariate across times t+1, t, and t-1 in the dataset.", call. = FALSE)
      } else {
        indcova2col <- which(names(data) == indcova[2])
      }
      
      if (length(indcova) == 3) {
        if (length(which(names(data) == indcova[3])) == 0 && indcova[3] != "none") {
          stop("Variable indcova must equal either NA or a vector of up to 3 exact names of an individual covariate across times t+1, t, and t-1 in the dataset.", call. = FALSE)
        } else {
          indcova1col <- which(names(data) == indcova[3])
        }
      } else {
        indcova[3] <- "none"
      }
    }
  } else {indcova <- c("none", "none", "none")}
  
  if (any(!is.na(indcovb))) {
    if (length(indcovb) > 3) {
      warning("Vector indcovb holds the names of an individual covariate across times t+1, t, and t-1. Will use only the first three elements.")
      indcovb <- indcovb[1:3]
    } else if (length(indcovb) == 1) {
      warning("Vector indcovb requires the names of an individual covariate across times t+1, t, and t-1. Since only 1 name was supplied, will assume that this name corresponds to status in time t+1 and not use the variable.")
      indcovb <- c("none", "none", "none")
    } else {
      if (length(which(names(data) == indcovb[2])) == 0 && indcovb[2] != "none") {
        stop("Variable indcovb must equal either NA or a vector of up to 3 exact names of an individual covariate across times t+1, t, and t-1 in the dataset.", call. = FALSE)
      } else {
        indcovb2col <- which(names(data) == indcovb[2])
      }
      
      if (length(indcovb) == 3) {
        if (length(which(names(data) == indcovb[3])) == 0 && indcovb[3] != "none") {
          stop("Variable indcovb must equal either NA or a vector of up to 3 exact names of an individual covariate across times t+1, t, and t-1 in the dataset.", call. = FALSE)
        } else {
          indcovb1col <- which(names(data) == indcovb[3])
        }
      } else {
        indcovb[3] <- "none"
      }
    }
  } else {indcovb <- c("none", "none", "none")}
  
  if (any(!is.na(indcovc))) {
    if (length(indcovc) > 3) {
      warning("Vector indcovc holds the names of an individual covariate across times t+1, t, and t-1. Will use only the first three elements.")
      indcovc <- indcovc[1:3]
    } else if (length(indcovc) == 1) {
      warning("Vector indcovc requires the names of an individual covariate across times t+1, t, and t-1. Since only 1 name was supplied, will assume that this name corresponds to status in time t+1 and not use the variable.")
      indcovc <- c("none", "none", "none")
    } else {
      if (length(which(names(data) == indcovc[2])) == 0 && indcovc[2] != "none") {
        stop("Variable indcovc must equal either NA or a vector of up to 3 exact names of an individual covariate across times t+1, t, and t-1 in the dataset.", call. = FALSE)
      } else {
        indcovc2col <- which(names(data) == indcovc[2])
      }
      
      if (length(indcovc) == 3) {
        if (length(which(names(data) == indcovc[3])) == 0 && indcovc[3] != "none") {
          stop("Variable indcovc must equal either NA or a vector of up to 3 exact names of an individual covariate across times t+1, t, and t-1 in the dataset.", call. = FALSE)
        } else {
          indcovc1col <- which(names(data) == indcovc[3])
        }
      } else {
        indcovc[3] <- "none"
      }
    }
  } else {indcovc <- c("none", "none", "none")}
  
  if (!is.na(indiv)) {
    if (length(which(names(data) == indiv)) == 0) {
      stop("Variable indiv must equal either NA or the exact name of the variable denoting individual identity in the dataset.", 
           call. = FALSE)
    } else {
      indivcol <- which(names(data) == indiv)
    }
  } else {indiv <- "none"}
  
  if (!is.na(patch)) {
    if (length(which(names(data) == patch)) == 0) {
      stop("Variable patch must equal either NA or the exact name of the variable denoting patch identity in the dataset.", 
           call. = FALSE)
    } else {
      patchcol <- which(names(data) == patch)
    }
  } else {patch <- "none"}
  
  if (!is.na(year)) {
    if (length(which(names(data) == year)) == 0) {
      stop("Variable year must equal either NA or the exact name of the variable denoting time t in the dataset.", 
           call. = FALSE)
    } else {
      yearcol <- which(names(data) == year)
    }
  } else {year <- "none"}
  
  # Here we test the dataset for appropriate stage names
  if (!is.na(juvestimate)) {
    if (!any(is.element(stage, names(data)))) {stop("Names of stage variables do not match dataset.", call. = FALSE)}
    
    stage3col <- which(names(data) == stage[1])
    stage2col <- which(names(data) == stage[2])
    if (length(stage) == 3) {stage1col <- which(names(data) == stage[3])} else {stage1col <- 0}
    
    if (!is.element(juvestimate, unique(data[,stage2col]))) {
      stop("Stage input in juvestimate not recognized in stage2 variable in dataset", call. = FALSE)
    }
  }
  
  if (!is.logical(juvsize)) {
    stop("Option juvsize must be set to either TRUE or FALSE.", call. = FALSE)
  }
  
  #Now we check whether the best-fit criterion is appropriate, and set another variable equal to what we need
  if (!is.element(tolower(bestfit), c("aicc", "aicc&k"))) {
    stop("Bestfit must equal either 'AICc' or 'AICc&k', with the latter as the default.", call. = FALSE)
  }
  used.criterion <- gsub("&k", "", bestfit) #This variable will be used once dredging is done to determine best-fit models
  
  #The next section creates the text lines needed for the main model calls, based on function input
  formulae <- stovokor(surv, obs, size, repst, fec, vitalrates, historical, suite, approach, sizedist,
                       fecdist, is.na(juvestimate), age, indcova, indcovb, indcovc, indiv, patch, year, 
                       patch.as.random, year.as.random, fectime, juvsize, size.zero, fec.zero)
  
  #Now we need to create the input datasets
  if (!all(is.na(censor))) {
    if (length(censor) == 1) {
      data$censor2 <- data[, which(names(data) == censor[1])]
    } else {
      data$censor3 <- data[, which(names(data) == censor[1])]
      data$censor2 <- data[, which(names(data) == censor[2])]
      if (length(censor) > 2) {data$censor1 <- data[, which(names(data) == censor[3])]}
    }
  } else {
    data$censor2 <- 1
  }
  
  data <- subset(data, censor2 == 1)
  
  if (!all(is.na(censor))) {
    if (length(censor) > 1) {
      data <- subset(data, censor3 == 1)
      if (length(censor) > 2) {
        data <- subset(data, censor1 == 1)
      }
    }
  }
  
  if (!is.na(juvestimate)) {
    juvindivs <- which(data[,stage2col] == juvestimate)
    adultindivs <- setdiff(c(1:length(data[,stage2col])), juvindivs)
    
    juvsurv.data <- subset(data, data[,stage2col] == juvestimate & data[,which(names(data) == surv[2])] == 1)
    juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
    juvsurv.trans <- dim(juvsurv.data)[1]
    
    if (suite == "full" | suite == "main" | suite == "size") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[2])]))) {
        warning("NAs in size variables may cause model selection to fail.")
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[3])]))) {
          warning("NAs in size variables may cause model selection to fail.")
        }
      }
    }
    
    if (suite == "full" | suite == "main" | suite == "rep") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[2])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.")
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[3])]))) {
          warning("NAs in reproductive status variables may cause model selection to fail.")
        }
      }
    }
    
    if (is.element(0, juvsurv.data$matstatus3)) {
      warning("Function modelsearch() assumes that all juveniles either die or transition to maturity within 1 year. Some individuals in this dataset appear to live longer as juveniles than assumptions allow.")
    }
    
    juvobs.data <- subset(juvsurv.data, juvsurv.data[, which(names(juvsurv.data) == surv[1])] == 1)
    juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
    juvobs.trans <- dim(juvobs.data)[1]
    
    if (formulae$full.obs.model != 1) {
      juvsize.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
    } else {
      juvsize.data <- juvobs.data
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
    }
    juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
    juvsize.trans <- dim(juvsize.data)[1]
    
    juvrepst.data <- juvsize.data
    juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
    juvrepst.trans <- dim(juvrepst.data)[1]
    
    if (dim(juvsurv.data)[1] < 100) {
      warning("Juvenile dataset is very small, and some models may fail given the size.")
    }
    
    data <- data[adultindivs,] #This line resets the main dataset to adults only
  }
  
  surv.data <- subset(data, data[,which(names(data) == surv[2])] == 1)
  surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
  surv.trans <- dim(surv.data)[1]
  
  if (suite == "full" | suite == "main" | suite == "size") {
    if (any(is.na(surv.data[, which(names(surv.data) == size[2])]))) {
      warning("NAs in size variables may cause model selection to fail.")
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == size[3])]))) {
        warning("NAs in size variables may cause model selection to fail.")
      }
    }
  }
  
  if (suite == "full" | suite == "main" | suite == "rep") {
    if (any(is.na(surv.data[, which(names(surv.data) == repst[2])]))) {
      warning("NAs in reproductive status variables may cause model selection to fail.")
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == repst[3])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.")
      }
    }
  }
  
  if (dim(surv.data)[1] < 100) {
    warning("Dataset is very small, and some models may fail given the size.")
  }
  
  if(any(!suppressWarnings(!is.na(as.numeric(as.character(surv.data[, which(names(surv.data) == size[1])])))))) {
    warning("Modelsearch(), flefko3(), flefko2(), and aflefko2() are made to work with numeric size variables. Use of categorical variables may result in errors and unexpected behavior.")
  }
  
  obs.data <- subset(surv.data, surv.data[, which(names(surv.data) == surv[1])] == 1)
  obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
  obs.trans <- dim(obs.data)[1]
  
  if (formulae$full.obs.model != 1) {
    size.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
  } else {
    size.data <- obs.data
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
  }
  size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
  size.trans <- dim(size.data)[1]
  
  repst.data <- size.data
  repst.ind <- length(unique(repst.data[, which(names(repst.data) == indiv)]))
  repst.trans <- dim(repst.data)[1]
  
  if (formulae$full.repst.model != 1) {
    fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[2])] == 1)
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
      fec.data <- fec.data[which(fec.data[, which(names(fec.data) == repst[2])] == 1),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
      fec.data <- fec.data[which(fec.data[, which(names(fec.data) == repst[1])] == 1),]
    }
  } else {
    fec.data <- surv.data
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
    }
  }
  fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
  fec.trans <- dim(fec.data)[1]
  
  #Now we check for exceptions to size in the dataset
  if (sizedist == "poisson") {
    
    if (any(size.data[, which(names(size.data) == size[1])] != round(size.data[, which(names(size.data) == size[1])]))) {
      stop("Size variables must be composed only of integers for the Poisson distribution to be used.")
    }
    
    if (!is.na(juvestimate)) {
      if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != round(juvsize.data[, which(names(juvsize.data) == size[1])]))) {
        stop("Size variables must be composed only of integers for the Poisson distribution to be used.")
      }
    }
    
  } else if (sizedist == "negbin") {
    
    if (any(size.data[, which(names(size.data) == size[1])] != round(size.data[, which(names(size.data) == size[1])]))) {
      stop("Size variables must be composed only of integers for the negative binomial distribution to be used.")
    }
    
    if (!is.na(juvestimate)) {
      if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != round(juvsize.data[, which(names(juvsize.data) == size[1])]))) {
        stop("Size variables must be composed only of integers for the negative binomial distribution to be used.")
      }
    }
  }
  
  if (fecdist == "poisson") {
    if (fectime == 2) {
      usedfec <- which(names(fec.data) == fec[2])
    } else if (fectime == 3) {
      usedfec <- which(names(fec.data) == fec[1])
    }
    
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers for the Poisson distribution to be used.")
    }
  } else if (fecdist == "negbin") {
    if (fectime == 2) {
      usedfec <- which(names(fec.data) == fec[2])
    } else if (fectime == 3) {
      usedfec <- which(names(fec.data) == fec[1])
    }
    
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers for the negative binomial distribution to be used.")
    }
  }
  
  #Now we run the modeling exercises
  if (formulae$full.surv.model == 1) {surv.global.model <- 1}
  if (formulae$full.obs.model == 1) {obs.global.model <- 1}
  if (formulae$full.size.model == 1) {size.global.model <- 1}
  if (formulae$full.repst.model == 1) {repst.global.model <- 1}
  if (formulae$full.fec.model == 1) {fec.global.model <- 1}
  
  if (formulae$juv.surv.model == 1) {juv.surv.global.model <- 1}
  if (formulae$juv.obs.model == 1) {juv.obs.global.model <- 1}
  if (formulae$juv.size.model == 1) {juv.size.global.model <- 1}
  if (formulae$juv.repst.model == 1) {juv.repst.global.model <- 1}
  
  if (formulae$juv.surv.model == 0) {juv.surv.global.model <- 0}
  if (formulae$juv.obs.model == 0) {juv.obs.global.model <- 0}
  if (formulae$juv.size.model == 0) {juv.size.global.model <- 0}
  if (formulae$juv.repst.model == 0) {juv.repst.global.model <- 0}
  
  surv.table <- NA
  obs.table <- NA
  size.table <- NA
  repst.table <- NA
  fec.table <- NA
  
  juvsurv.table <- NA
  juvobs.table <- NA
  juvsize.table <- NA
  juvrepst.table <- NA
  
  surv.bf <- NA
  obs.bf <- NA
  size.bf <- NA
  repst.bf <- NA
  fec.bf <- NA
  
  juvsurv.bf <- NA
  juvobs.bf <- NA
  juvsize.bf <- NA
  juvrepst.bf <- NA
  
  #A few more corrections to the model structure, used in running the global models
  correction.sz1sz2 <- gsub("sz1", size[3], " + sz1:sz2", fixed = TRUE)
  correction.sz1sz2 <- gsub("sz2", size[2], correction.sz1sz2, fixed = TRUE)
  correction.fl1fl2 <- gsub("fl1", repst[3], " + fl1:fl2", fixed = TRUE)
  correction.fl1fl2 <- gsub("fl2", repst[2], correction.fl1fl2, fixed = TRUE)
  correction.sz1fl1 <- gsub("sz1", size[3], " + sz1:fl1", fixed = TRUE)
  correction.sz1fl1 <- gsub("fl1", repst[3], correction.sz1fl1, fixed = TRUE)
  correction.sz2fl2 <- gsub("sz2", size[2], " + sz2:fl2", fixed = TRUE)
  correction.sz2fl2 <- gsub("fl2", repst[2], correction.sz2fl2, fixed = TRUE)
  correction.sz1fl2 <- gsub("sz1", size[3], " + sz1:fl2", fixed = TRUE)
  correction.sz1fl2 <- gsub("fl2", repst[2], correction.sz1fl2, fixed = TRUE)
  correction.sz2fl1 <- gsub("sz2", size[2], " + sz2:fl1", fixed = TRUE)
  correction.sz2fl1 <- gsub("fl1", repst[3], correction.sz2fl1, fixed = TRUE)
  
  correction.sz1sz2a <- gsub("sz1", size[3], " + sz2:sz1", fixed = TRUE)
  correction.sz1sz2a <- gsub("sz2", size[2], correction.sz1sz2a, fixed = TRUE)
  correction.fl1fl2a <- gsub("fl1", repst[3], " + fl2:fl1", fixed = TRUE)
  correction.fl1fl2a <- gsub("fl2", repst[2], correction.fl1fl2a, fixed = TRUE)
  correction.sz1fl1a <- gsub("sz1", size[3], " + fl1:sz1", fixed = TRUE)
  correction.sz1fl1a <- gsub("fl1", repst[3], correction.sz1fl1a, fixed = TRUE)
  correction.sz2fl2a <- gsub("sz2", size[2], " + fl2:sz2", fixed = TRUE)
  correction.sz2fl2a <- gsub("fl2", repst[2], correction.sz2fl2a, fixed = TRUE)
  correction.sz1fl2a <- gsub("sz1", size[3], " + fl2:sz1", fixed = TRUE)
  correction.sz1fl2a <- gsub("fl2", repst[2], correction.sz1fl2a, fixed = TRUE)
  correction.sz2fl1a <- gsub("sz2", size[2], " + fl1:sz2", fixed = TRUE)
  correction.sz2fl1a <- gsub("fl1", repst[3], correction.sz2fl1a, fixed = TRUE)
  
  correction.indiv <- gsub("individ", indiv, " + (1 | individ)", fixed = TRUE)
  if(approach == "mixed" & year.as.random) {
    correction.year <- gsub("yr", year, " + (1 | yr)", fixed = TRUE)
  } else {
    correction.year <- gsub("yr", year, " + yr", fixed = TRUE)
  }
  if (approach == "mixed" & patch.as.random) {
    correction.patch <- gsub("patch", patch, " + (1 | patch)", fixed = TRUE)
  } else {
    correction.patch <- gsub("patch", patch, " + patch", fixed = TRUE)
  }
  
  #Here we run the global models
  if (approach == "mixed") {
    
    if (formulae$full.surv.model != 1) {
      
      if (is.element(0, surv.data$alive3) & is.element(1, surv.data$alive3)) {
        
        if (!quiet) {writeLines("\nDeveloping global model of survival probability...\n"); }
        surv.global.model <- try(lme4::glmer(formula = stats::as.formula(formulae$full.surv.model), data = surv.data, family = "binomial"),
                                 silent = TRUE)
        
        if (any(class(surv.global.model) == "try-error")) {
          nox.surv.model <- formulae$full.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.surv.model <- gsub(correction.sz1sz2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz1sz2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.surv.model <- gsub(correction.fl1fl2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.fl1fl2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.surv.model <- gsub(correction.sz1fl1, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz1fl1a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.surv.model <- gsub(correction.sz2fl2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz2fl2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.surv.model <- gsub(correction.sz1fl2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz1fl2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.surv.model <- gsub(correction.sz2fl1, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz2fl1a, "", nox.surv.model, fixed = TRUE)
          }
          
          if (nox.surv.model != formulae$full.surv.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            surv.global.model <- try(lme4::glmer(formula = nox.surv.model, data = surv.data, family = "binomial"),
                                     silent = TRUE)
          }
          
          nopat.surv.model <- nox.surv.model
          if (!is.na(correction.patch)) {
            nopat.surv.model <- gsub(correction.patch, "", nopat.surv.model, fixed = TRUE)
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            
            if (nox.surv.model != nopat.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              surv.global.model <- try(lme4::glmer(formula = nopat.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noyr.surv.model <- nopat.surv.model
          if (!is.na(correction.year)) {
            noyr.surv.model <- gsub(correction.year, "", noyr.surv.model, fixed = TRUE)
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            
            if (noyr.surv.model != nopat.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              surv.global.model <- try(lme4::glmer(formula = noyr.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noind.surv.model <- noyr.surv.model
          if (!is.na(correction.indiv)) {
            noind.surv.model <- gsub(correction.indiv, "", noind.surv.model, fixed = TRUE)
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            
            if (noind.surv.model != noyr.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              surv.global.model <- try(lme4::glmer(formula = noind.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for survival.")}
            
            surv.global.model <- 1
            surv.ind <- 0
            surv.trans <- 0
          }
        }
      } else if (!is.element(0, surv.data$alive3)) {
        if (!quiet) {writeLines("\nSurvival response is constant so will not model it.")}
        formulae$full.surv.model <- 1
        surv.global.model <- 1
        surv.ind <- 0
        surv.trans <- 0
      } else {
        if (!quiet) {writeLines("\nSurvival response is constant so will not model it.")}
        formulae$full.surv.model <- 1
        surv.global.model <- 0
        surv.ind <- 0
        surv.trans <- 0
      }
    } else {
      surv.global.model <- 1
      surv.ind <- 0
      surv.trans <- 0
    }
    
    if (formulae$full.obs.model != 1) {
      if (is.element(0, obs.data$obsstatus3) & is.element(1, obs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of observation probability...\n"); }
        obs.global.model <- try(lme4::glmer(formula = formulae$full.obs.model, data = obs.data, family = "binomial"),
                                silent = TRUE)
        
        if (any(class(obs.global.model) == "try-error")) {
          nox.obs.model <- formulae$full.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.obs.model <- gsub(correction.sz1sz2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz1sz2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.obs.model <- gsub(correction.fl1fl2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.fl1fl2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.obs.model <- gsub(correction.sz1fl1, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz1fl1a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.obs.model <- gsub(correction.sz2fl2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz2fl2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.obs.model <- gsub(correction.sz1fl2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz1fl2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.obs.model <- gsub(correction.sz2fl1, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz2fl1a, "", nox.obs.model, fixed = TRUE)
          }
          
          if (nox.obs.model != formulae$full.obs.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            obs.global.model <- try(lme4::glmer(formula = nox.obs.model, data = obs.data, family = "binomial"),
                                    silent = TRUE)
          }
          
          nopat.obs.model <- nox.obs.model
          if (!is.na(correction.patch)) {
            nopat.obs.model <- gsub(correction.patch, "", nopat.obs.model, fixed = TRUE)
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            
            if (nox.obs.model != nopat.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              obs.global.model <- try(lme4::glmer(formula = nopat.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noyr.obs.model <- nopat.obs.model
          if (!is.na(correction.year)) {
            noyr.obs.model <- gsub(correction.year, "", noyr.obs.model, fixed = TRUE)
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            
            if (noyr.obs.model != nopat.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              obs.global.model <- try(lme4::glmer(formula = noyr.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noind.obs.model <- noyr.obs.model
          if (!is.na(correction.indiv)) {
            noind.obs.model <- gsub(correction.indiv, "", noyr.obs.model, fixed = TRUE)
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            
            if (noind.obs.model != noyr.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              obs.global.model <- try(lme4::glmer(formula = noind.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            if (!quiet) {writeLines("Could not properly estimate a global model for observation status.")}
            
            obs.global.model <- 1
            obs.ind <- 0
            obs.trans <- 0
          }
        }
      } else if (!is.element(0, obs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nObservation response is constant so will not model it.")}
        formulae$full.obs.model <- 1
        obs.global.model <- 1
        obs.ind <- 0
        obs.trans <- 0
      } else {
        if (!quiet) {writeLines("\nObservation response is constant so will not model it.")}
        formulae$full.obs.model <- 1
        obs.global.model <- 0
        obs.ind <- 0
        obs.trans <- 0
      }
    } else {
      obs.global.model <- 1
      obs.ind <- 0
      obs.trans <- 0
    }
    
    if (formulae$full.size.model != 1) {
      if (sizedist == "gaussian") {
        if (!quiet) {writeLines("\nDeveloping global model of size (Gaussian)...\n");}
        size.global.model <- try(lme4::lmer(formula = formulae$full.size.model, data = size.data),
                                 silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(lme4::lmer(formula = nox.size.model, data = size.data), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(lme4::lmer(formula = nopat.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(lme4::lmer(formula = noyr.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(lme4::lmer(formula = noind.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "poisson" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (Poisson)...\n");}
        size.global.model <- try(lme4::glmer(formula = formulae$full.size.model, data = size.data, family = "poisson"),
                                 silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(lme4::glmer(formula = nox.size.model, data = size.data, family = "poisson"),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(lme4::glmer(formula = nopat.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(lme4::glmer(formula = noyr.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(lme4::glmer(formula = noind.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (negative binomial)...\n");}
        size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$full.size.model), data = size.data,
                                                  ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.size.model), data = size.data,
                                                      ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.size.model), data = size.data,
                                                        ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.size.model), data = size.data,
                                                        ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.size.model), data = size.data,
                                                        ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      }  else if (sizedist == "poisson" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (Poisson)...\n");}
        size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$full.size.model), data = size.data,
                                                  ziformula=~., family = "poisson"), silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.size.model), data = size.data,
                                                      ziformula=~., family = "poisson"), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.size.model), data = size.data,
                                                        ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.size.model), data = size.data,
                                                        ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.size.model), data = size.data,
                                                        ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (negative binomial)...\n");}
        size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$full.size.model), data = size.data,
                                                  ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.size.model), data = size.data,
                                                      ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.size.model), data = size.data,
                                                        ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.size.model), data = size.data,
                                                        ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.size.model), data = size.data,
                                                        ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      }
    } else {
      size.global.model <- 1
      size.ind <- 0
      size.trans <- 0
    }
    
    if (formulae$full.repst.model != 1) {
      if (is.element(0, repst.data$repstatus3) & is.element(1, repst.data$repstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of the probability of reproduction...\n"); }
        repst.global.model <- try(lme4::glmer(formula = formulae$full.repst.model, data = repst.data, family = "binomial"),
                                  silent = TRUE)
        
        if (any(class(repst.global.model) == "try-error")) {
          nox.repst.model <- formulae$full.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.repst.model <- gsub(correction.sz1sz2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz1sz2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.repst.model <- gsub(correction.fl1fl2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.fl1fl2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.repst.model <- gsub(correction.sz1fl1, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz1fl1a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.repst.model <- gsub(correction.sz2fl2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz2fl2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.repst.model <- gsub(correction.sz1fl2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz1fl2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.repst.model <- gsub(correction.sz2fl1, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz2fl1a, "", nox.repst.model, fixed = TRUE)
          }
          
          if (nox.repst.model != formulae$full.repst.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            repst.global.model <- try(lme4::glmer(formula = nox.repst.model, data = repst.data, family = "binomial"),
                                      silent = TRUE)
          }
          
          nopat.repst.model <- nox.repst.model
          if (!is.na(correction.patch)) {
            nopat.repst.model <- gsub(correction.patch, "", nopat.repst.model, fixed = TRUE)
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            
            if (nox.repst.model != nopat.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              repst.global.model <- try(lme4::glmer(formula = nopat.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noyr.repst.model <- nopat.repst.model
          if (!is.na(correction.year)) {
            noyr.repst.model <- gsub(correction.year, "", noyr.repst.model, fixed = TRUE)
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            
            if (noyr.repst.model != nopat.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              repst.global.model <- try(lme4::glmer(formula = noyr.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noind.repst.model <- noyr.repst.model
          if (!is.na(correction.indiv)) {
            noind.repst.model <- gsub(correction.indiv, "", noind.repst.model, fixed = TRUE)
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            
            if (noind.repst.model != noyr.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              repst.global.model <- try(lme4::glmer(formula = noind.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for reproduction status.")}
            repst.global.model <- 1
            repst.ind <- 0
            repst.trans <- 0
          }
        }
      } else if (!is.element(0, repst.data$repstatus3)) {
        if (!quiet) {writeLines("\nReproductive status response is constant so will not model it.")}
        formulae$full.repst.model <- 1
        repst.global.model <- 1
        repst.ind <- 0
        repst.trans <- 0
      } else {
        if (!quiet) {writeLines("\nReproductive status response is constant so will not model it.")}
        formulae$full.repst.model <- 1
        repst.global.model <- 0
        repst.ind <- 0
        repst.trans <- 0
      }
    } else {
      repst.global.model <- 1
      repst.ind <- 0
      repst.trans <- 0
    }
    
    if (formulae$full.fec.model != 1) {
      if (fecdist == "gaussian") {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (Gaussian)...\n");}
        fec.global.model <- try(lme4::lmer(formula = formulae$full.fec.model, data = fec.data), silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(lme4::lmer(formula = nox.fec.model, data = fec.data),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(lme4::lmer(formula = nopat.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.patch)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(lme4::lmer(formula = noyr.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(lme4::lmer(formula = noind.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "poisson" & !fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (Poisson)...\n");}
        fec.global.model <- try(lme4::glmer(formula = formulae$full.fec.model, data = fec.data, family = "poisson"),
                                silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(lme4::glmer(formula = nox.fec.model, data = fec.data, family = "poisson"),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(lme4::glmer(formula = nopat.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(lme4::glmer(formula = noyr.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(lme4::glmer(formula = noind.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "negbin" & !fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (negative binomial)...\n");}
        fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$full.fec.model), data = fec.data,
                                                 ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.fec.model), data = fec.data,
                                                     ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.fec.model), data = fec.data,
                                                       ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.fec.model), data = fec.data,
                                                       ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.fec.model), data = fec.data,
                                                       ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      }  else if (fecdist == "poisson" & fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (Poisson)...\n");}
        fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$full.fec.model), data = fec.data,
                                                 ziformula=~., family = "poisson"), silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.fec.model), data = fec.data,
                                                     ziformula=~., family = "poisson"), silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.fec.model), data = fec.data,
                                                       ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.fec.model), data = fec.data,
                                                       ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.fec.model), data = fec.data,
                                                       ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "negbin" & fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (negative binomial)...\n");}
        fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$full.fec.model), data = fec.data,
                                                 ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.fec.model), data = fec.data,
                                                     ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.fec.model), data = fec.data,
                                                       ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.fec.model), data = fec.data,
                                                       ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.fec.model), data = fec.data,
                                                       ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      }
    } else {
      fec.global.model <- 1
      fec.ind <- 0
      fec.trans <- 0
    }
    
    if (formulae$juv.surv.model != 1 & formulae$juv.surv.model != 0) {
      if (is.element(0, juvsurv.data$alive3) & is.element(1, juvsurv.data$alive3)) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile survival probability...\n"); }
        juv.surv.global.model <- try(lme4::glmer(formula = formulae$juv.surv.model, data = juvsurv.data, family = "binomial"),
                                     silent = TRUE)
        
        if (any(class(juv.surv.global.model) == "try-error")) {
          nox.juv.surv.model <- formulae$juv.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.surv.model <- gsub(correction.sz1sz2, "", nox.juv.surv.model, fixed = TRUE)
            nox.juv.surv.model <- gsub(correction.sz1sz2a, "", nox.juv.surv.model, fixed = TRUE)
          }
          
          if (nox.juv.surv.model != formulae$juv.surv.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.surv.global.model <- try(lme4::glmer(formula = nox.juv.surv.model, data = juvsurv.data, 
                                                     family = "binomial"), silent = TRUE)
          }
          
          nopat.juv.surv.model <- nox.juv.surv.model
          if (!is.na(correction.patch)) {
            nopat.juv.surv.model <- gsub(correction.patch, "", nopat.juv.surv.model, fixed = TRUE)
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            
            if (nox.juv.surv.model != nopat.juv.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.surv.global.model <- try(lme4::glmer(formula = nopat.juv.surv.model, data = juvsurv.data, 
                                                       family = "binomial"), silent = TRUE)
            }
          }
          
          noyr.juv.surv.model <- nopat.juv.surv.model
          if (!is.na(correction.year)) {
            noyr.juv.surv.model <- gsub(correction.year, "", noyr.juv.surv.model, fixed = TRUE)
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            
            if (noyr.juv.surv.model != nopat.juv.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.surv.global.model <- try(lme4::glmer(formula = noyr.juv.surv.model, data = juvsurv.data, 
                                                       family = "binomial"), silent = TRUE)
            }
          }
          
          noind.juv.surv.model <- noyr.juv.surv.model
          if (!is.na(correction.indiv)) {
            noind.juv.surv.model <- gsub(correction.indiv, "", noind.juv.surv.model, fixed = TRUE)
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            
            if (noind.juv.surv.model != noyr.juv.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.surv.global.model <- try(lme4::glmer(formula = noind.juv.surv.model, data = juvsurv.data, 
                                                       family = "binomial"), silent = TRUE)
            }
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile survival.")}
            juv.surv.global.model <- 1
            juvsurv.ind <- 0
            juvsurv.trans <- 0
          }
        }
      } else if (!is.element(0, juvsurv.data$alive3)) {
        if (!quiet) {writeLines("\nJuvenile survival response is constant so will not model it.")}
        formulae$juv.surv.model <- 1
        juv.surv.global.model <- 1
        juvsurv.ind <- 0
        juvsurv.trans <- 0
      } else {
        if (!quiet) {writeLines("\nJuvenile survival response is constant so will not model it.")}
        formulae$juv.surv.model <- 1
        juv.surv.global.model <- 0
        juvsurv.ind <- 0
        juvsurv.trans <- 0
      }
    } else {
      juv.surv.global.model <- 1
      juvsurv.ind <- 0
      juvsurv.trans <- 0
    }
    
    if (formulae$juv.obs.model != 1 & formulae$juv.obs.model != 0) {
      if (is.element(0, juvobs.data$obsstatus3) & is.element(1, juvobs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile observation probability...\n"); }
        juv.obs.global.model <- try(lme4::glmer(formula = formulae$juv.obs.model, data = juvobs.data, 
                                                family = "binomial"), silent = TRUE)
        
        if (any(class(juv.obs.global.model) == "try-error")) {
          nox.juv.obs.model <- formulae$juv.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.obs.model <- gsub(correction.sz1sz2, "", nox.juv.obs.model, fixed = TRUE)
            nox.juv.obs.model <- gsub(correction.sz1sz2a, "", nox.juv.obs.model, fixed = TRUE)
          }
          
          if (nox.juv.obs.model != formulae$juv.obs.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.obs.global.model <- try(lme4::glmer(formula = nox.juv.obs.model, data = juvobs.data, 
                                                    family = "binomial"), silent = TRUE)
          }
          
          nopat.juv.obs.model <- nox.juv.obs.model
          if (!is.na(correction.patch)) {
            nopat.juv.obs.model <- gsub(correction.patch, "", nopat.juv.obs.model, fixed = TRUE)
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            
            if (nox.juv.obs.model != nopat.juv.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.obs.global.model <- try(lme4::glmer(formula = nopat.juv.obs.model, data = juvobs.data, 
                                                      family = "binomial"), silent = TRUE)
            }
          }
          
          noyr.juv.obs.model <- nopat.juv.obs.model
          if (!is.na(correction.year)) {
            noyr.juv.obs.model <- gsub(correction.year, "", noyr.juv.obs.model, fixed = TRUE)
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            
            if (noyr.juv.obs.model != nopat.juv.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.obs.global.model <- try(lme4::glmer(formula = noyr.juv.obs.model, data = juvobs.data, 
                                                      family = "binomial"), silent = TRUE)
            }
          }
          
          noind.juv.obs.model <- noyr.juv.obs.model
          if (!is.na(correction.indiv)) {
            noind.juv.obs.model <- gsub(correction.indiv, "", noind.juv.obs.model, fixed = TRUE)
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            
            if (noind.juv.obs.model != noyr.juv.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.obs.global.model <- try(lme4::glmer(formula = noind.juv.obs.model, data = juvobs.data, 
                                                      family = "binomial"), silent = TRUE)
            }
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile observation status.")}
            juv.obs.global.model <- 1
            juvobs.ind <- 0
            juvobs.trans <- 0
          }
        }
      } else if (!is.element(0, juvobs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nJuvenile observation response is constant so will not model it.")}
        formulae$juv.obs.model <- 1
        juv.obs.global.model <- 1
        juvobs.ind <- 0
        juvobs.trans <- 0
      } else {
        if (!quiet) {writeLines("\nJuvenile observation response is constant so will not model it.")}
        formulae$juv.obs.model <- 1
        juv.obs.global.model <- 0
        juvobs.ind <- 0
        juvobs.trans <- 0
      }
    } else {
      juv.obs.global.model <- 1
      juvobs.ind <- 0
      juvobs.trans <- 0
    }
    
    if (formulae$juv.size.model != 1 & formulae$juv.size.model != 0) {
      if (sizedist == "gaussian") {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile size (Gaussian)...\n");}
        juv.size.global.model <- try(lme4::lmer(formula = formulae$juv.size.model, data = juvsize.data),
                                     silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(lme4::lmer(formula = nox.juv.size.model, data = juvsize.data), 
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(lme4::lmer(formula = nopat.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(lme4::lmer(formula = noyr.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(lme4::lmer(formula = noind.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      } else if (sizedist == "poisson" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile size (Poisson)...\n");}
        juv.size.global.model <- try(lme4::glmer(formula = formulae$juv.size.model, data = juvsize.data, 
                                                 family = "poisson"), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(lme4::glmer(formula = nox.juv.size.model, data = juvsize.data, 
                                                     family = "poisson"), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("Global model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(lme4::glmer(formula = nopat.juv.size.model, data = juvsize.data, 
                                                       family = "poisson"), silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(lme4::glmer(formula = noyr.juv.size.model, data = juvsize.data, 
                                                       family = "poisson"), silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(lme4::glmer(formula = noind.juv.size.model, data = juvsize.data, 
                                                       family = "poisson"), silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile size (negative binomial)...\n");}
        juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$juv.size.model), data = juvsize.data, 
                                                      ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.juv.size.model), data = juvsize.data, ziformula=~0, 
                                                          family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.juv.size.model), 
                                                            data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.juv.size.model), 
                                                            data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.juv.size.model), 
                                                            data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                           silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      }  else if (sizedist == "poisson" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile size (Poisson)...\n");}
        juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$juv.size.model), data = juvsize.data, 
                                                      ziformula=~., family = "poisson"), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.juv.size.model), data = juvsize.data, ziformula=~., 
                                                          family = "poisson"), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.juv.size.model), data = juvsize.data, 
                                                            ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.juv.size.model), data = juvsize.data, 
                                                            ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.juv.size.model), data = juvsize.data, 
                                                            ziformula=~., family = "poisson"), silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile size (negative binomial)...\n");}
        juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(formulae$juv.size.model), data = juvsize.data, 
                                                      ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nox.juv.size.model), data = juvsize.data, ziformula=~., 
                                                          family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(nopat.juv.size.model), data = juvsize.data, 
                                                            ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noyr.juv.size.model), data = juvsize.data, 
                                                            ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(noind.juv.size.model), data = juvsize.data, 
                                                            ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      }
    } else {
      juv.size.global.model <- 1
      juvsize.ind <- 0
      juvsize.trans <- 0
    }
    
    if (formulae$juv.repst.model != 1 & formulae$juv.repst.model != 0) {
      if (is.element(0, juvrepst.data$repstatus3) & is.element(1, juvrepst.data$repstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile reproduction probability...\n"); }
        juv.repst.global.model <- try(lme4::glmer(formula = formulae$juv.repst.model, data = juvrepst.data, 
                                                  family = "binomial"), silent = TRUE)
        
        if (any(class(juv.repst.global.model) == "try-error")) {
          nox.juv.repst.model <- formulae$juv.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.repst.model <- gsub(correction.sz1sz2, "", nox.juv.repst.model, fixed = TRUE)
            nox.juv.repst.model <- gsub(correction.sz1sz2a, "", nox.juv.repst.model, fixed = TRUE)
          }
          
          if (nox.juv.repst.model != formulae$juv.repst.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.repst.global.model <- try(lme4::glmer(formula = nox.juv.repst.model, data = juvrepst.data, 
                                                      family = "binomial"), silent = TRUE)
          }
          
          nopat.juv.repst.model <- nox.juv.repst.model
          if (!is.na(correction.patch)) {
            nopat.juv.repst.model <- gsub(correction.patch, "", nopat.juv.repst.model, fixed = TRUE)
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            
            if (nox.juv.repst.model != nopat.juv.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.repst.global.model <- try(lme4::glmer(formula = nopat.juv.repst.model, data = juvrepst.data,
                                                        family = "binomial"), silent = TRUE)
            }
          }
          
          noyr.juv.repst.model <- nopat.juv.repst.model
          if (!is.na(correction.year)) {
            noyr.juv.repst.model <- gsub(correction.year, "", noyr.juv.repst.model, fixed = TRUE)
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            
            if (noyr.juv.repst.model != nopat.juv.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.repst.global.model <- try(lme4::glmer(formula = noyr.juv.repst.model, data = juvrepst.data,
                                                        family = "binomial"), silent = TRUE)
            }
          }
          
          noind.juv.repst.model <- noyr.juv.repst.model
          if (!is.na(correction.indiv)) {
            noind.juv.repst.model <- gsub(correction.indiv, "", noind.juv.repst.model, fixed = TRUE)
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            
            if (noind.juv.repst.model != noyr.juv.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.repst.global.model <- try(lme4::glmer(formula = noind.juv.repst.model, data = juvrepst.data,
                                                        family = "binomial"), silent = TRUE)
            }
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile reproduction status.")}
            juv.repst.global.model <- 1
            juvrepst.ind <- 0
            juvrepst.trans <- 0
          }
        }
      } else if (!is.element(0, juvrepst.data$repstatus3)) {
        if (!quiet) {writeLines("\nJuvenile reproduction response is constant so will not model it.")}
        formulae$juv.repst.model <- 1
        juv.repst.global.model <- 1
        juvrepst.ind <- 0
        juvrepst.trans <- 0
      } else {
        if (!quiet) {writeLines("\nJuvenile reproduction response is constant so will not model it.")}
        formulae$juv.repst.model <- 1
        juv.repst.global.model <- 0
        juvrepst.ind <- 0
        juvrepst.trans <- 0
      }
    } else {
      juvrepst.ind <- 0
      juvrepst.trans <- 0
    }
    
    if (!quiet) {writeLines("\nAll global models developed.\n")}
    
  } else if (approach == "glm") {
    
    if (formulae$full.surv.model != 1) {
      
      if (is.element(0, surv.data$alive3) & is.element(1, surv.data$alive3)) {
        
        if (!quiet) {writeLines("\nDeveloping global model of survival probability...\n"); }
        surv.global.model <- try(glm(formula = formulae$full.surv.model, data = surv.data, family = "binomial"),
                                 silent = TRUE)
        
        if (any(class(surv.global.model) == "try-error")) {
          nox.surv.model <- formulae$full.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.surv.model <- gsub(correction.sz1sz2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz1sz2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.surv.model <- gsub(correction.fl1fl2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.fl1fl2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.surv.model <- gsub(correction.sz1fl1, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz1fl1a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.surv.model <- gsub(correction.sz2fl2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz2fl2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.surv.model <- gsub(correction.sz1fl2, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz1fl2a, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.surv.model <- gsub(correction.sz2fl1, "", nox.surv.model, fixed = TRUE)
            nox.surv.model <- gsub(correction.sz2fl1a, "", nox.surv.model, fixed = TRUE)
          }
          
          if (nox.surv.model != formulae$full.surv.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            surv.global.model <- try(glm(formula = nox.surv.model, data = surv.data, family = "binomial"),
                                     silent = TRUE)
          }
          
          nopat.surv.model <- nox.surv.model
          if (!is.na(correction.patch)) {
            nopat.surv.model <- gsub(correction.patch, "", nopat.surv.model, fixed = TRUE)
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            
            if (nox.surv.model != nopat.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              surv.global.model <- try(glm(formula = nopat.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noyr.surv.model <- nopat.surv.model
          if (!is.na(correction.year)) {
            noyr.surv.model <- gsub(correction.year, "", noyr.surv.model, fixed = TRUE)
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            
            if (noyr.surv.model != nopat.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              surv.global.model <- try(glm(formula = noyr.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noind.surv.model <- noyr.surv.model
          if (!is.na(correction.indiv)) {
            noind.surv.model <- gsub(correction.indiv, "", noind.surv.model, fixed = TRUE)
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            
            if (noind.surv.model != noyr.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              surv.global.model <- try(glm(formula = noind.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          if (any(class(surv.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for survival.")}
            
            surv.global.model <- 1
            surv.ind <- 0
            surv.trans <- 0
          }
        }
      } else if (!is.element(0, surv.data$alive3)) {
        if (!quiet) {writeLines("\nSurvival response is constant so will not model it.")}
        formulae$full.surv.model <- 1
        surv.global.model <- 1
        surv.ind <- 0
        surv.trans <- 0
      } else {
        if (!quiet) {writeLines("\nSurvival response is constant so will not model it.")}
        formulae$full.surv.model <- 1
        surv.global.model <- 0
        surv.ind <- 0
        surv.trans <- 0
      }
    } else {
      surv.global.model <- 1
      surv.ind <- 0
      surv.trans <- 0
    }
    
    if (formulae$full.obs.model != 1) {
      if (is.element(0, obs.data$obsstatus3) & is.element(1, obs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of observation probability...\n"); }
        obs.global.model <- try(glm(formula = formulae$full.obs.model, data = obs.data, family = "binomial"),
                                silent = TRUE)
        
        if (any(class(obs.global.model) == "try-error")) {
          nox.obs.model <- formulae$full.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.obs.model <- gsub(correction.sz1sz2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz1sz2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.obs.model <- gsub(correction.fl1fl2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.fl1fl2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.obs.model <- gsub(correction.sz1fl1, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz1fl1a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.obs.model <- gsub(correction.sz2fl2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz2fl2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.obs.model <- gsub(correction.sz1fl2, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz1fl2a, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.obs.model <- gsub(correction.sz2fl1, "", nox.obs.model, fixed = TRUE)
            nox.obs.model <- gsub(correction.sz2fl1a, "", nox.obs.model, fixed = TRUE)
          }
          
          if (nox.obs.model != formulae$full.obs.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            obs.global.model <- try(glm(formula = nox.obs.model, data = obs.data, family = "binomial"),
                                    silent = TRUE)
          }
          
          nopat.obs.model <- nox.obs.model
          if (!is.na(correction.patch)) {
            nopat.obs.model <- gsub(correction.patch, "", nopat.obs.model, fixed = TRUE)
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            
            if (nox.obs.model != nopat.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              obs.global.model <- try(glm(formula = nopat.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noyr.obs.model <- nopat.obs.model
          if (!is.na(correction.year)) {
            noyr.obs.model <- gsub(correction.year, "", noyr.obs.model, fixed = TRUE)
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            
            if (noyr.obs.model != nopat.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              obs.global.model <- try(glm(formula = noyr.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noind.obs.model <- noyr.obs.model
          if (!is.na(correction.indiv)) {
            noind.obs.model <- gsub(correction.indiv, "", noind.obs.model, fixed = TRUE)
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            
            if (noind.obs.model != noyr.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              obs.global.model <- try(glm(formula = noind.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          if (any(class(obs.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for observation status.")}
            
            obs.global.model <- 1
            obs.ind <- 0
            obs.trans <- 0
          }
        }
      } else if (!is.element(0, obs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nObservation response is constant so will not model it.")}
        formulae$full.obs.model <- 1
        obs.global.model <- 1
        obs.ind <- 0
        obs.trans <- 0
      } else {
        if (!quiet) {writeLines("\nObservation response is constant so will not model it.")}
        formulae$full.obs.model <- 1
        obs.global.model <- 0
        obs.ind <- 0
        obs.trans <- 0
      }
    } else {
      obs.global.model <- 1
      obs.ind <- 0
      obs.trans <- 0
    }
    
    if (formulae$full.size.model != 1) {
      if (sizedist == "gaussian") {
        if (!quiet) {writeLines("\nDeveloping global model of size (Gaussian)...\n");}
        size.global.model <- try(lm(formula = formulae$full.size.model, data = size.data), silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(lm(formula = nox.size.model, data = size.data), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(lm(formula = nopat.size.model, data = size.data), silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(lm(formula = noyr.size.model, data = size.data), silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(lm(formula = noind.size.model, data = size.data), silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "poisson" & !size.zero) {
        
        if (!quiet) {writeLines("\nDeveloping global model of size (Poisson)...\n");}
        size.global.model <- try(glm(formula = formulae$full.size.model, data = size.data, family = "poisson"),
                                 silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(glm(formula = nox.size.model, data = size.data, family = "poisson"),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(glm(formula = nopat.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(glm(formula = noyr.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(glm(formula = noind.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (negative binomial)...\n");}
        size.global.model <- try(MASS::glm.nb(formula = formulae$full.size.model, data = size.data), silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(MASS::glm.nb(formula = nox.size.model, data = size.data),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(MASS::glm.nb(formula = nopat.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(MASS::glm.nb(formula = noyr.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(MASS::glm.nb(formula = noind.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "poisson" & size.zero) {
        
        if (!quiet) {writeLines("\nDeveloping global model of size (Poisson)...\n");}
        size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(formulae$full.size.model), data = size.data, dist = "poisson"),
                                 silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nox.size.model), data = size.data, dist = "poisson"),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nopat.size.model), data = size.data, dist = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noyr.size.model), data = size.data, dist = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noind.size.model), data = size.data, dist = "poisson"),
                                       silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (negative binomial)...\n");}
        size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(formulae$full.size.model), 
                                                data = size.data, dist = "negbin"), silent = TRUE)
        
        if (any(class(size.global.model) == "try-error")) {
          nox.size.model <- formulae$full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1sz2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.fl1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl1a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz1fl2a, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
            nox.size.model <- gsub(correction.sz2fl1a, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != formulae$full.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nox.size.model), data = size.data, dist = "negbin"),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (nox.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nopat.size.model), data = size.data, dist = "negbin"),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noyr.size.model != nopat.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noyr.size.model), data = size.data, dist = "negbin"),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (any(class(size.global.model) == "try-error")) {
            
            if (noind.size.model != noyr.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noind.size.model), data = size.data, dist = "negbin"),
                                       silent = TRUE)
            }
          }
          
          if (any(class(size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for size.")}
            
            size.global.model <- 1
            size.ind <- 0
            size.trans <- 0
          }
        }
      }
    } else {
      size.global.model <- 1
      size.ind <- 0
      size.trans <- 0
    }
    
    
    if (formulae$full.repst.model != 1) {
      if (is.element(0, repst.data$repstatus3) & is.element(1, repst.data$repstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of the probability of reproduction...\n"); }
        repst.global.model <- try(glm(formula = formulae$full.repst.model, data = repst.data, family = "binomial"),
                                  silent = TRUE)
        
        if (any(any(class(repst.global.model) == "try-error"))) {
          nox.repst.model <- formulae$full.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.repst.model <- gsub(correction.sz1sz2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz1sz2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.repst.model <- gsub(correction.fl1fl2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.fl1fl2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.repst.model <- gsub(correction.sz1fl1, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz1fl1a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.repst.model <- gsub(correction.sz2fl2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz2fl2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.repst.model <- gsub(correction.sz1fl2, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz1fl2a, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.repst.model <- gsub(correction.sz2fl1, "", nox.repst.model, fixed = TRUE)
            nox.repst.model <- gsub(correction.sz2fl1a, "", nox.repst.model, fixed = TRUE)
          }
          
          if (nox.repst.model != formulae$full.repst.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            repst.global.model <- try(glm(formula = nox.repst.model, data = repst.data, family = "binomial"),
                                      silent = TRUE)
          }
          
          nopat.repst.model <- nox.repst.model
          if (!is.na(correction.patch)) {
            nopat.repst.model <- gsub(correction.patch, "", nopat.repst.model, fixed = TRUE)
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            
            if (nox.repst.model != nopat.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              repst.global.model <- try(glm(formula = nopat.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noyr.repst.model <- nopat.repst.model
          if (!is.na(correction.year)) {
            noyr.repst.model <- gsub(correction.year, "", noyr.repst.model, fixed = TRUE)
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            
            if (noyr.repst.model != nopat.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              repst.global.model <- try(glm(formula = noyr.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noind.repst.model <- noyr.repst.model
          if (!is.na(correction.indiv)) {
            noind.repst.model <- gsub(correction.indiv, "", noind.repst.model, fixed = TRUE)
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            
            if (noind.repst.model != noyr.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              repst.global.model <- try(glm(formula = noind.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          if (any(class(repst.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for reproduction status.")}
            
            repst.global.model <- 1
            repst.ind <- 0
            repst.trans <- 0
          }
        }
      } else if (!is.element(0, repst.data$repstatus3)) {
        if (!quiet) {writeLines("\nReproductive status response is constant so will not model it.")}
        formulae$full.repst.model <- 1
        repst.global.model <- 1
        repst.ind <- 0
        repst.trans <- 0
      } else {
        if (!quiet) {writeLines("\nReproductive status response is constant so will not model it.")}
        formulae$full.repst.model <- 1
        repst.global.model <- 0
        repst.ind <- 0
        repst.trans <- 0
      }
    } else {
      repst.global.model <- 1
      repst.ind <- 0
      repst.trans <- 0
    }
    
    if (formulae$full.fec.model != 1) {
      if (fecdist == "gaussian") {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (Gaussian)...\n");}
        fec.global.model <- try(lm(formula = formulae$full.fec.model, data = fec.data), silent = TRUE)
        
        if (any(any(class(fec.global.model) == "try-error"))) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(lm(formula = nox.fec.model, data = fec.data), silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(lm(formula = nopat.fec.model, data = fec.data), silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(lm(formula = noyr.fec.model, data = fec.data), silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(lm(formula = noind.fec.model, data = fec.data), silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "poisson" & !fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (Poisson)...\n");}
        fec.global.model <- try(glm(formula = formulae$full.fec.model, data = fec.data, family = "poisson"),
                                silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(glm(formula = nox.fec.model, data = fec.data, family = "poisson"),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(glm(formula = nopat.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(glm(formula = noyr.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(glm(formula = noind.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "negbin" & !fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (negative binomial)...\n");}
        fec.global.model <- try(MASS::glm.nb(formula = formulae$full.fec.model, data = fec.data), silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(MASS::glm.nb(formula = nox.fec.model, data = fec.data),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(MASS::glm.nb(formula = nopat.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(MASS::glm.nb(formula = noyr.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(MASS::glm.nb(formula = noind.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "poisson" & fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (Poisson)...\n");}
        fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(formulae$full.fec.model), data = fec.data, dist = "poisson"),
                                silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nox.fec.model), data = fec.data, dist = "poisson"),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nopat.fec.model), data = fec.data, dist = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noyr.fec.model), data = fec.data, dist = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noind.fec.model), data = fec.data, dist = "poisson"),
                                      silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      } else if (fecdist == "negbin" & fec.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of fecundity (negative binomial)...\n");}
        fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(formulae$full.fec.model), data = fec.data, dist = "negbin"), 
                                silent = TRUE)
        
        if (any(class(fec.global.model) == "try-error")) {
          nox.fec.model <- formulae$full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1sz2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.fl1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl1a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz1fl2a, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
            nox.fec.model <- gsub(correction.sz2fl1a, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != formulae$full.fec.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nox.fec.model), data = fec.data, dist = "negbin"),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (nox.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nopat.fec.model), data = fec.data, dist = "negbin"),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noyr.fec.model != nopat.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noyr.fec.model), data = fec.data, dist = "negbin"),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            
            if (noind.fec.model != noyr.fec.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              fec.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noind.fec.model), data = fec.data, dist = "negbin"),
                                      silent = TRUE)
            }
          }
          
          if (any(class(fec.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for fecundity.")}
            fec.global.model <- 1
            fec.ind <- 0
            fec.trans <- 0
          }
        }
      }
    } else {
      fec.global.model <- 1
      fec.ind <- 0
      fec.trans <- 0
    }
    
    if (formulae$juv.surv.model != 1 & formulae$juv.surv.model != 0) {
      if (is.element(0, juvsurv.data$alive3) & is.element(1, juvsurv.data$alive3)) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile survival probability...\n"); }
        juv.surv.global.model <- try(glm(formula = formulae$juv.surv.model, data = juvsurv.data, family = "binomial"),
                                     silent = TRUE)
        
        if (any(class(juv.surv.global.model) == "try-error")) {
          nox.juv.surv.model <- formulae$juv.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.surv.model <- gsub(correction.sz1sz2, "", nox.juv.surv.model, fixed = TRUE)
            nox.juv.surv.model <- gsub(correction.sz1sz2a, "", nox.juv.surv.model, fixed = TRUE)
          }
          
          if (nox.juv.surv.model != formulae$juv.surv.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.surv.global.model <- try(glm(formula = nox.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                         silent = TRUE)
          }
          
          nopat.juv.surv.model <- nox.juv.surv.model
          if (!is.na(correction.patch)) {
            nopat.juv.surv.model <- gsub(correction.patch, "", nopat.juv.surv.model, fixed = TRUE)
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            
            if (nox.juv.surv.model != nopat.juv.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.surv.global.model <- try(glm(formula = nopat.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.surv.model <- nopat.juv.surv.model
          if (!is.na(correction.year)) {
            noyr.juv.surv.model <- gsub(correction.year, "", noyr.juv.surv.model, fixed = TRUE)
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            
            if (noyr.juv.surv.model != nopat.juv.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.surv.global.model <- try(glm(formula = noyr.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                           silent = TRUE)
            }
          }
          
          noind.juv.surv.model <- noyr.juv.surv.model
          if (!is.na(correction.indiv)) {
            noind.juv.surv.model <- gsub(correction.indiv, "", noind.juv.surv.model, fixed = TRUE)
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            
            if (noind.juv.surv.model != noyr.juv.surv.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.surv.global.model <- try(glm(formula = noind.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                           silent = TRUE)
            }
          }
          
          if (any(class(juv.surv.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile survival.")}
            juv.surv.global.model <- 1
            juvsurv.ind <- 0
            juvsurv.trans <- 0
          }
        }
      } else if (!is.element(0, juvsurv.data$alive3)) {
        if (!quiet) {writeLines("\nJuvenile survival response is constant so will not model it.")}
        formulae$juv.surv.model <- 1
        juv.surv.global.model <- 1
        juvsurv.ind <- 0
        juvsurv.trans <- 0
      } else {
        if (!quiet) {writeLines("\nJuvenile survival response is constant so will not model it.")}
        formulae$juv.surv.model <- 1
        juv.surv.global.model <- 0
        juvsurv.ind <- 0
        juvsurv.trans <- 0
      }
    } else {
      juv.surv.global.model <- 1
      juvsurv.ind <- 0
      juvsurv.trans <- 0
    }
    
    if (formulae$juv.obs.model != 1 & formulae$juv.obs.model != 0) {
      if (is.element(0, juvobs.data$obsstatus3) & is.element(1, juvobs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile observation probability...\n"); }
        juv.obs.global.model <- try(glm(formula = formulae$juv.obs.model, data = juvobs.data, family = "binomial"),
                                    silent = TRUE)
        
        if (any(class(juv.obs.global.model) == "try-error")) {
          nox.juv.obs.model <- formulae$juv.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.obs.model <- gsub(correction.sz1sz2, "", nox.juv.obs.model, fixed = TRUE)
            nox.juv.obs.model <- gsub(correction.sz1sz2a, "", nox.juv.obs.model, fixed = TRUE)
          }
          
          if (nox.juv.obs.model != formulae$juv.obs.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.obs.global.model <- try(glm(formula = nox.juv.obs.model, data = juvobs.data, family = "binomial"),
                                        silent = TRUE)
          }
          
          nopat.juv.obs.model <- nox.juv.obs.model
          if (!is.na(correction.patch)) {
            nopat.juv.obs.model <- gsub(correction.patch, "", nopat.juv.obs.model, fixed = TRUE)
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            
            if (nox.juv.obs.model != nopat.juv.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.obs.global.model <- try(glm(formula = nopat.juv.obs.model, data = juvobs.data, family = "binomial"),
                                          silent = TRUE)
            }
          }
          
          noyr.juv.obs.model <- nopat.juv.obs.model
          if (!is.na(correction.year)) {
            noyr.juv.obs.model <- gsub(correction.year, "", noyr.juv.obs.model, fixed = TRUE)
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            
            if (noyr.juv.obs.model != nopat.juv.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.obs.global.model <- try(glm(formula = noyr.juv.obs.model, data = juvobs.data, family = "binomial"),
                                          silent = TRUE)
            }
          }
          
          noind.juv.obs.model <- noyr.juv.obs.model
          if (!is.na(correction.indiv)) {
            noind.juv.obs.model <- gsub(correction.indiv, "", noind.juv.obs.model, fixed = TRUE)
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            
            if (noind.juv.obs.model != noyr.juv.obs.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.obs.global.model <- try(glm(formula = noind.juv.obs.model, data = juvobs.data, family = "binomial"),
                                          silent = TRUE)
            }
          }
          
          if (any(class(juv.obs.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile observation status.")}
            juv.obs.global.model <- 1
            juvobs.ind <- 0
            juvobs.trans <- 0
          }
        }
      } else if (!is.element(0, juvobs.data$obsstatus3)) {
        if (!quiet) {writeLines("\nJuvenile observation response is constant so will not model it.")}
        formulae$juv.obs.model <- 1
        juv.obs.global.model <- 1
        juvobs.ind <- 0
        juvobs.trans <- 0
      } else {
        if (!quiet) {writeLines("\nJuvenile observation response is constant so will not model it.")}
        formulae$juv.obs.model <- 1
        juv.obs.global.model <- 0
        juvobs.ind <- 0
        juvobs.trans <- 0
      }
    } else {
      juv.obs.global.model <- 1
      juvobs.ind <- 0
      juvobs.trans <- 0
    }
    
    if (formulae$juv.size.model != 1 & formulae$juv.size.model != 0) {
      if (sizedist == "gaussian") {
        if (!quiet) {writeLines("\nDeveloping global model of size (Gaussian)...\n");}
        juv.size.global.model <- try(lm(formula = formulae$juv.size.model, data = juvsize.data), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(lm(formula = nox.juv.size.model, data = juvsize.data),
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(lm(formula = nopat.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(lm(formula = noyr.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(lm(formula = noind.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      } else if (sizedist == "poisson" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (Poisson)...\n");}
        juv.size.global.model <- try(glm(formula = formulae$juv.size.model, data = juvsize.data, family = "poisson"),
                                     silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(glm(formula = nox.juv.size.model, data = juvsize.data, family = "poisson"),
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(glm(formula = nopat.juv.size.model, data = juvsize.data, family = "poisson"),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(glm(formula = noyr.juv.size.model, data = juvsize.data, family = "poisson"),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(glm(formula = noind.juv.size.model, data = juvsize.data, family = "poisson"),
                                           silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & !size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (negative binomial)...\n");}
        juv.size.global.model <- try(MASS::glm.nb(formula = formulae$juv.size.model, data = juvsize.data),
                                     silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(MASS::glm.nb(formula = nox.juv.size.model, data = juvsize.data),
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(MASS::glm.nb(formula = nopat.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(MASS::glm.nb(formula = noyr.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(MASS::glm.nb(formula = noind.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      }  else if (sizedist == "poisson" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (Poisson)...\n");}
        juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(formulae$juv.size.model), 
                                                    data = juvsize.data, dist = "poisson"), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nox.juv.size.model), 
                                                        data = juvsize.data, dist = "poisson"), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nopat.juv.size.model), 
                                                          data = juvsize.data, dist = "poisson"), silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noyr.juv.size.model), 
                                                          data = juvsize.data, dist = "poisson"), silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noind.juv.size.model), 
                                                          data = juvsize.data, dist = "poisson"), silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      } else if (sizedist == "negbin" & size.zero) {
        if (!quiet) {writeLines("\nDeveloping global model of size (negative binomial)...\n");}
        juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(formulae$juv.size.model), 
                                                    data = juvsize.data, dist = "negbin"), silent = TRUE)
        
        if (any(class(juv.size.global.model) == "try-error")) {
          nox.juv.size.model <- formulae$juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
            nox.juv.size.model <- gsub(correction.sz1sz2a, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != formulae$juv.size.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nox.juv.size.model), 
                                                        data = juvsize.data, dist = "negbin"), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(nopat.juv.size.model), 
                                                          data = juvsize.data, dist = "negbin"), silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noyr.juv.size.model), 
                                                          data = juvsize.data, dist = "negbin"), silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.size.global.model <- try(pscl::zeroinfl(formula = stats::as.formula(noind.juv.size.model), 
                                                          data = juvsize.data, dist = "negbin"), silent = TRUE)
            }
          }
          
          if (any(class(juv.size.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile size.")}
            juv.size.global.model <- 1
            juvsize.ind <- 0
            juvsize.trans <- 0
          }
        }
      }
    } else {
      juv.size.global.model <- 1
      juvsize.ind <- 0
      juvsize.trans <- 0
    }
    
    if (formulae$juv.repst.model != 1 & formulae$juv.repst.model != 0) {
      if (is.element(0, juvrepst.data$repstatus3) & is.element(1, juvrepst.data$repstatus3)) {
        if (!quiet) {writeLines("\nDeveloping global model of juvenile reproduction probability...\n"); }
        juv.repst.global.model <- try(glm(formula = formulae$juv.repst.model, data = juvrepst.data, family = "binomial"),
                                      silent = TRUE)
        
        if (any(class(juv.repst.global.model) == "try-error")) {
          nox.juv.repst.model <- formulae$juv.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.repst.model <- gsub(correction.sz1sz2, "", nox.juv.repst.model, fixed = TRUE)
            nox.juv.repst.model <- gsub(correction.sz1sz2a, "", nox.juv.repst.model, fixed = TRUE)
          }
          
          if (nox.juv.repst.model != formulae$juv.repst.model) {
            if (!quiet) {writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")}
            juv.repst.global.model <- try(glm(formula = nox.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                          silent = TRUE)
          }
          
          nopat.juv.repst.model <- nox.juv.repst.model
          if (!is.na(correction.patch)) {
            nopat.juv.repst.model <- gsub(correction.patch, "", nopat.juv.repst.model, fixed = TRUE)
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            
            if (nox.juv.repst.model != nopat.juv.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")}
              juv.repst.global.model <- try(glm(formula = nopat.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                            silent = TRUE)
            }
          }
          
          noyr.juv.repst.model <- nopat.juv.repst.model
          if (!is.na(correction.year)) {
            noyr.juv.repst.model <- gsub(correction.year, "", noyr.juv.repst.model, fixed = TRUE)
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            
            if (noyr.juv.repst.model != nopat.juv.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")}
              juv.repst.global.model <- try(glm(formula = noyr.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                            silent = TRUE)
            }
          }
          
          noind.juv.repst.model <- noyr.juv.repst.model
          if (!is.na(correction.indiv)) {
            noind.juv.repst.model <- gsub(correction.indiv, "", noind.juv.repst.model, fixed = TRUE)
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            
            if (noind.juv.repst.model != noyr.juv.repst.model) {
              if (!quiet) {writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")}
              juv.repst.global.model <- try(glm(formula = noind.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                            silent = TRUE)
            }
          }
          
          if (any(class(juv.repst.global.model) == "try-error")) {
            if (!quiet) {writeLines("\nCould not properly estimate a global model for juvenile observation status.")}
            juv.repst.global.model <- 1
            juvrepst.ind <- 0
            juvrepst.trans <- 0
          }
        }
      } else if (!is.element(0, juvrepst.data$repstatus3)) {
        if (!quiet) {writeLines("\nJuvenile observation response is constant so will not model it.")}
        formulae$juv.repst.model <- 1
        juv.repst.global.model <- 1
        juvrepst.ind <- 0
        juvrepst.trans <- 0
      } else {
        if (!quiet) {writeLines("\nJuvenile observation response is constant so will not model it.")}
        formulae$juv.repst.model <- 1
        juv.repst.global.model <- 0
        juvrepst.ind <- 0
        juvrepst.trans <- 0
      }
    } else {
      juvrepst.ind <- 0
      juvrepst.trans <- 0
    }
    
    if (!quiet) {writeLines("All global models developed.\n")}
  }
  
  #This is the section where we dredge the models
  if (suite != "cons" & global.only == FALSE) {
    options(na.action = "na.fail");
    
    if (formulae$full.surv.model != 1) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of survival probability...\n")
        surv.table <- try(MuMIn::dredge(surv.global.model, rank = used.criterion), silent = TRUE)
      } else {
        surv.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(surv.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(surv.table) == "try-error")) {warning("Dredge of survival probability failed.")}
    }
    
    if (formulae$full.obs.model != 1) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of observation probability...\n")
        obs.table <- try(MuMIn::dredge(obs.global.model, rank = used.criterion), silent = TRUE)
      } else {
        obs.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(obs.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(obs.table) == "try-error")) {warning("Dredge of observation probability failed.")}
    }
    
    if (formulae$full.size.model != 1) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of size...\n")
        size.table <- try(MuMIn::dredge(size.global.model, rank = used.criterion), silent = TRUE)
      } else {
        size.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(size.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(size.table) == "try-error")) {warning("Dredge of size response failed.")}
    }
    
    if (formulae$full.repst.model != 1) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of reproduction probability...\n")
        repst.table <- try(MuMIn::dredge(repst.global.model, rank = used.criterion), silent = TRUE)
      } else {
        repst.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(repst.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(repst.table) == "try-error")) {warning("Dredge of reproductive status probability failed.")}
    }
    
    if (formulae$full.fec.model != 1) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of fecundity...\n")
        fec.table <- try(MuMIn::dredge(fec.global.model, rank = used.criterion), silent = TRUE)
      } else {
        fec.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(fec.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(fec.table) == "try-error")) {warning("Dredge of fecundity response failed.")}
    }
    
    if (formulae$juv.surv.model != 1 & formulae$juv.surv.model != 0) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of juvenile survival probability...\n")
        juvsurv.table <- try(MuMIn::dredge(juv.surv.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvsurv.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(juv.surv.global.model, rank = used.criterion), 
                                                               silent = TRUE)))
      }
      if (any(class(juvsurv.table) == "try-error")) {warning("Dredge of juvenilesurvival probability failed.")}
    }
    
    if (formulae$juv.obs.model != 1 & formulae$juv.obs.model != 0) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of juvenile observation probability...\n")
        juvobs.table <- try(MuMIn::dredge(juv.obs.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvobs.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(juv.obs.global.model, rank = used.criterion), 
                                                              silent = TRUE)))
      }
      if (any(class(juvobs.table) == "try-error")) {warning("Dredge of juvenile observation probability failed.")}
    }
    
    if (formulae$juv.size.model != 1 & formulae$juv.size.model != 0) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of juvenile size...\n")
        juvsize.table <- try(MuMIn::dredge(juv.size.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvsize.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(juv.size.global.model, rank = used.criterion), 
                                                               silent = TRUE)))
      }
      if (any(class(juvsize.table) == "try-error")) {warning("Dredge of juvenile size response failed.")}
    }
    
    if (formulae$juv.repst.model != 1 & formulae$juv.repst.model != 0) {
      options(na.action = "na.fail")
      if (!quiet) {
        writeLines("\nCommencing dredge of juvenile reproduction status...\n")
        juvrepst.table <- try(MuMIn::dredge(juv.repst.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvrepst.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(juv.repst.global.model, rank = used.criterion), 
                                                                silent = TRUE)))
      }
      if (any(class(juvrepst.table) == "try-error")) {warning("Dredge of juvenile reproduction status failed.")}
    }
    
    if (!quiet) {writeLines("\nFinished dredging all vital rates.\n")}
  }
  
  if (stringr::str_detect(bestfit, "&k")) {
    if (any(class(surv.table) == "model.selection")) {
      relevant.surv.models <- which(surv.table$delta <= 2)
      df.surv.models <- which(surv.table$df[relevant.surv.models] == min(surv.table$df[relevant.surv.models]))
      if (length(df.surv.models) > 1) {
        df.surv.models <- min(df.surv.models)
      }
      if (quiet) {
        surv.bf <- suppressWarnings(suppressMessages(eval(getCall(surv.table, min(df.surv.models)))))
      } else {
        surv.bf <- eval(getCall(surv.table, min(df.surv.models)))
      }
    } else {
      surv.bf <- surv.global.model
    }
    
    if (any(class(obs.table) == "model.selection")) {
      relevant.obs.models <- which(obs.table$delta <= 2)
      df.obs.models <- which(obs.table$df[relevant.obs.models] == min(obs.table$df[relevant.obs.models]))
      if (length(df.obs.models) > 1) {
        df.obs.models <- min(df.obs.models)
      }
      if (quiet) {
        obs.bf <- suppressWarnings(suppressMessages(eval(getCall(obs.table, min(df.obs.models)))))
      } else {
        obs.bf <- eval(getCall(obs.table, min(df.obs.models)))
      }
    } else {
      obs.bf <- obs.global.model
    }
    
    if (any(class(size.table) == "model.selection")) {
      relevant.size.models <- which(size.table$delta <= 2)
      df.size.models <- which(size.table$df[relevant.size.models] == min(size.table$df[relevant.size.models]))
      if (length(df.size.models) > 1) {
        df.size.models <- min(df.size.models)
      }
      if (quiet) {
        size.bf <- suppressWarnings(suppressMessages(eval(getCall(size.table, min(df.size.models)))))
      } else {
        size.bf <- eval(getCall(size.table, min(df.size.models)))
      }
    } else {
      size.bf <- size.global.model
    }
    
    if (any(class(repst.table) == "model.selection")) {
      relevant.repst.models <- which(repst.table$delta <= 2)
      df.repst.models <- which(repst.table$df[relevant.repst.models] == min(repst.table$df[relevant.repst.models]))
      if (length(df.repst.models) > 1) {
        df.repst.models <- min(df.repst.models)
      }
      if (quiet) {
        repst.bf <- suppressWarnings(suppressMessages(eval(getCall(repst.table, min(df.repst.models)))))
      } else {
        repst.bf <- eval(getCall(repst.table, min(df.repst.models)))
      }
    } else {
      repst.bf <- repst.global.model
    }
    
    if (any(class(fec.table) == "model.selection")) {
      relevant.fec.models <- which(fec.table$delta <= 2)
      df.fec.models <- which(fec.table$df[relevant.fec.models] == min(fec.table$df[relevant.fec.models]))
      if (length(df.fec.models) > 1) {
        df.fec.models <- min(df.fec.models)
      }
      if (quiet) {
        fec.bf <- suppressWarnings(suppressMessages(eval(getCall(fec.table, min(df.fec.models)))))
      } else {
        fec.bf <- eval(getCall(fec.table, min(df.fec.models)))
      }
    } else {
      fec.bf <- fec.global.model
    }
    
    if (any(class(juvsurv.table) == "model.selection")) {
      relevant.juv.surv.models <- which(juvsurv.table$delta <= 2)
      df.juv.surv.models <- which(juvsurv.table$df[relevant.juv.surv.models] == min(juvsurv.table$df[relevant.juv.surv.models]))
      if (length(df.juv.surv.models) > 1) {
        df.juv.surv.models <- min(df.juv.surv.models)
      }
      if (quiet) {
        juvsurv.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsurv.table, min(df.juv.surv.models)))))
      } else {
        juvsurv.bf <- eval(getCall(juvsurv.table, min(df.juv.surv.models)))
      }
    } else {
      juvsurv.bf <- juv.surv.global.model
    }
    
    if (any(class(juvobs.table) == "model.selection")) {
      relevant.juv.obs.models <- which(juvobs.table$delta <= 2)
      df.juv.obs.models <- which(juvobs.table$df[relevant.juv.obs.models] == min(juvobs.table$df[relevant.juv.obs.models]))
      if (length(df.juv.obs.models) > 1) {
        df.juv.obs.models <- min(df.juv.obs.models)
      }
      if (quiet) {
        juvobs.bf <- suppressWarnings(suppressMessages(eval(getCall(juvobs.table, min(df.juv.obs.models)))))
      } else {
        juvobs.bf <- eval(getCall(juvobs.table, min(df.juv.obs.models)))
      }
    } else {
      juvobs.bf <- juv.obs.global.model
    }
    
    if (any(class(juvsize.table) == "model.selection")) {
      relevant.juv.size.models <- which(juvsize.table$delta <= 2)
      df.juv.size.models <- which(juvsize.table$df[relevant.juv.size.models] == min(juvsize.table$df[relevant.juv.size.models]))
      if (length(df.juv.size.models) > 1) {
        df.juv.size.models <- min(df.juv.size.models)
      }
      if (quiet) {
        juvsize.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsize.table, min(df.juv.size.models)))))
      } else {
        juvsize.bf <- eval(getCall(juvsize.table, min(df.juv.size.models)))
      }
    } else {
      juvsize.bf <- juv.size.global.model
    }
    
    if (any(class(juvrepst.table) == "model.selection")) {
      relevant.juv.repst.models <- which(juvrepst.table$delta <= 2)
      df.juv.repst.models <- which(juvrepst.table$df[relevant.juv.repst.models] == min(juvrepst.table$df[relevant.juv.repst.models]))
      if (length(df.juv.repst.models) > 1) {
        df.juv.repst.models <- min(df.juv.repst.models)
      }
      if (quiet) {
        juvrepst.bf <- suppressWarnings(suppressMessages(eval(getCall(juvrepst.table, min(df.juv.repst.models)))))
      } else {
        juvrepst.bf <- eval(getCall(juvrepst.table, min(df.juv.repst.models)))
      }
    } else {
      juvrepst.bf <- juv.repst.global.model
    }
  } else if (!stringr::str_detect(bestfit, "&k")) {
    if (any(class(surv.table) == "model.selection")) {
      if (quiet) {
        surv.bf <- suppressWarnings(suppressMessages(eval(getCall(surv.table, 1))))
      } else {
        surv.bf <- eval(getCall(surv.table, 1))
      }
    } else {
      surv.bf <- surv.global.model
    }
    
    if (any(class(obs.table) == "model.selection")) {
      if (quiet) {
        obs.bf <- suppressWarnings(suppressMessages(eval(getCall(obs.table, 1))))
      } else {
        obs.bf <- eval(getCall(obs.table, 1))
      }
    } else {
      obs.bf <- obs.global.model
    }
    
    if (any(class(size.table) == "model.selection")) {
      if (quiet) {
        size.bf <- suppressWarnings(suppressMessages(eval(getCall(size.table, 1))))
      } else {
        size.bf <- eval(getCall(size.table, 1))
      }
    } else {
      size.bf <- size.global.model
    }
    
    if (any(class(repst.table) == "model.selection")) {
      if (quiet) {
        repst.bf <- suppressWarnings(suppressMessages(eval(getCall(repst.table, 1))))
      } else {
        repst.bf <- eval(getCall(repst.table, 1))
      }
    } else {
      repst.bf <- repst.global.model
    }
    
    if (any(class(fec.table) == "model.selection")) {
      if (quiet) {
        fec.bf <- suppressWarnings(suppressMessages(eval(getCall(fec.table, 1))))
      } else {
        fec.bf <- eval(getCall(fec.table, 1))
      }
    } else {
      fec.bf <- fec.global.model
    }
    
    if (any(class(juvsurv.table) == "model.selection")) {
      if (quiet) {
        juvsurv.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsurv.table, 1))))
      } else {
        juvsurv.bf <- eval(getCall(juvsurv.table, 1))
      }
    } else {
      juvsurv.bf <- juv.surv.global.model
    }
    
    if (any(class(juvobs.table) == "model.selection")) {
      if (quiet) {
        juvobs.bf <- suppressWarnings(suppressMessages(eval(getCall(juvobs.table, 1))))
      } else {
        juvobs.bf <- eval(getCall(juvobs.table, 1))
      }
    } else {
      juvobs.bf <- juv.obs.global.model
    }
    
    if (any(class(juvsize.table) == "model.selection")) {
      if (quiet) {
        juvsize.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsize.table, 1))))
      } else {
        juvsize.bf <- eval(getCall(juvsize.table, 1))
      }
    } else {
      juvsize.bf <- juv.size.global.model
    }
    
    if (any(class(juvrepst.table) == "model.selection")) {
      if (quiet) {
        juvrepst.bf <- suppressWarnings(suppressMessages(eval(getCall(juvrepst.table, 1))))
      } else {
        juvrepst.bf <- eval(getCall(juvrepst.table, 1))
      }
    } else {
      juvrepst.bf <- juv.repst.global.model
    }
  }
  
  if (!quiet & global.only == FALSE) {writeLines("\nFinished selecting best-fit models.\n")}
  
  qcoutput <- cbind.data.frame(c("survival", "observation", "size", "reproduction", "fecundity", 
                                 "juvenile_survival", "juvnile_observation", "juvenile_size", 
                                 "juvenile_reproduction"), c(surv.ind, obs.ind, size.ind, repst.ind, 
                                                             fec.ind, juvsurv.ind, juvobs.ind, juvsize.ind, juvrepst.ind), 
                               c(surv.trans, obs.trans, size.trans, repst.trans, fec.trans, 
                                 juvsurv.trans, juvobs.trans, juvsize.trans, juvrepst.trans))
  names(qcoutput) <- c("vital_rate", "individuals", "transitions")
  
  if (show.model.tables == FALSE) {
    surv.table <- NA
    obs.table <- NA
    size.table <- NA
    repst.table <- NA
    fec.table <- NA
    
    juvsurv.table <- NA
    juvobs.table <- NA
    juvsize.table <- NA
    juvrepst.table <- NA
  }
  
  if (global.only) {bestfit <- "global model only"}
  
  #Now we develop the final output, creating a new S3 class to do it
  full.output <- list(survival_model = surv.bf, observation_model = obs.bf, size_model = size.bf, 
                      repstatus_model = repst.bf, fecundity_model = fec.bf, 
                      juv_survival_model = juvsurv.bf, juv_observation_model = juvobs.bf, 
                      juv_size_model = juvsize.bf, juv_reproduction_model = juvrepst.bf,
                      survival_table = surv.table, observation_table = obs.table, 
                      size_table = size.table, repstatus_table = repst.table,
                      fecundity_table = fec.table, juv_survival_table = juvsurv.table, 
                      juv_observation_table = juvobs.table, juv_size_table = juvsize.table,
                      juv_reproduction_table = juvrepst.table, paramnames = formulae$paramnames, 
                      criterion = bestfit, qc = qcoutput)
  class(full.output) <- "lefkoMod"
  
  return(full.output)
}

#' Summary of Class "lefkoMod"
#' 
#' A function to sumarize the viewable output for an R object of class \code{lefkoMod}.
#' This function shows the best-fit models, summarizes the numbers of models in the
#' model tables, shows the criterion used to determine the best-fit models, and
#' provides some basic quality control information.
#' 
#' @param object An R object of class \code{lefkoMod} resulting from \code{\link{modelsearch}()}.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the best-fit models for all vital
#' rates, with constants of 0 or 1 used for unestimated models. This is followed by
#' a summary of the number of models tested per vital rate, and a table showing the
#' names of the parameters used to model vital rates and represent tested factors.
#' At the end is a section describing the number of individuals and individual
#' transitions used to estimate each vital rate best-fit model.
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
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "lme4", suite = "main",
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl", bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE, quiet = TRUE)
#' 
#' summary(lathmodelsln2)
#' }
#' 
#' @export
summary.lefkoMod <- function(object, ...) {
  
  modelsuite <- object

  totalmodels <- length(which(c(is.numeric(modelsuite$survival_model), is.numeric(modelsuite$observation_model), 
                                is.numeric(modelsuite$size_model), is.numeric(modelsuite$repstatus_model), 
                                is.numeric(modelsuite$fecundity_model), is.numeric(modelsuite$juv_survival_model),
                                is.numeric(modelsuite$juv_observation_model), is.numeric(modelsuite$juv_size_model),
                                is.numeric(modelsuite$juv_reproduction_model)) == FALSE))
  
  writeLines(paste0("This LefkoMod object includes ", totalmodels, " linear models."))
  writeLines(paste0("Best-fit model criterion used: ", modelsuite$criterion))
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("Survival model:")
  print(modelsuite$survival_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nObservation model:")
  print(modelsuite$observation_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nSize model:")
  print(modelsuite$size_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nReproductive status model:")
  print(modelsuite$repstatus_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nFecundity model:")
  print(modelsuite$fecundity_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("Juvenile survival model:")
  print(modelsuite$juv_survival_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nJuvenile observation model:")
  print(modelsuite$juv_observation_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nJuvenile size model:")
  print(modelsuite$juv_size_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nJuvenile reproduction model:")
  print(modelsuite$juv_reproduction_model)
  
  writeLines("\n\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  if (!is.logical(modelsuite$survival_model) & is.element("model.selection", class(modelsuite$survival_table))) {
    writeLines(paste0("\nNumber of models in survival table:", dim(modelsuite$survival_table)[1]))
  } else if (!is.logical(modelsuite$survival_model)) {
    writeLines("\nNumber of models in survival table: 1")
  } else {
    writeLines("\nSurvival table not estimated")
  }
  
  if (!is.logical(modelsuite$observation_model) & is.element("model.selection", class(modelsuite$observation_table))) {
    writeLines(paste0("\nNumber of models in observation table:", dim(modelsuite$observation_table)[1]))
  } else if (!is.logical(modelsuite$observation_model)) {
    writeLines("\nNumber of models in observation table: 1")
  } else {
    writeLines("\nObservation table not estimated")
  }
  
  if (!is.logical(modelsuite$size_model) & is.element("model.selection", class(modelsuite$size_table))) {
    writeLines(paste0("\nNumber of models in size table:", dim(modelsuite$size_table)[1]))
  } else if (!is.logical(modelsuite$size_model)) {
    writeLines("\nNumber of models in size table: 1")
  } else {
    writeLines("\nSize table not estimated")
  }
  
  if (!is.logical(modelsuite$repstatus_model) & is.element("model.selection", class(modelsuite$repstatus_table))) {
    writeLines(paste0("\nNumber of models in reproduction status table:", dim(modelsuite$repstatus_table)[1]))
  } else if (!is.logical(modelsuite$repstatus_model)) {
    writeLines("\nNumber of models in reproduction status table: 1")
  } else {
    writeLines("\nReproduction status table not estimated")
  }
  
  if (!is.logical(modelsuite$fecundity_model) & is.element("model.selection", class(modelsuite$fecundity_table))) {
    writeLines(paste0("\nNumber of models in fecundity table:", dim(modelsuite$fecundity_table)[1]))
  } else if (!is.logical(modelsuite$fecundity_model)) {
    writeLines("\nNumber of models in fecundity table: 1")
  } else {
    writeLines("\nFecundity table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_survival_model) & is.element("model.selection", class(modelsuite$juv_survival_table))) {
    writeLines(paste0("\nNumber of models in juvenile survival table:", dim(modelsuite$juv_survival_table)[1]))
  } else if (!is.logical(modelsuite$juv_survival_model)) {
    writeLines("\nNumber of models in juvenile survival table: 1")
  } else {
    writeLines("\nJuvenile survival table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_observation_model) & is.element("model.selection", class(modelsuite$juv_observation_table))) {
    writeLines(paste0("\nNumber of models in juvenile observation table:", dim(modelsuite$juv_observation_table)[1]))
  } else if (!is.logical(modelsuite$juv_observation_model)) {
    writeLines("\nNumber of models in juvenile observation table: 1")
  } else {
    writeLines("\nJuvenile observation table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_size_model) & is.element("model.selection", class(modelsuite$juv_size_table))) {
    writeLines(paste0("\nNumber of models in juvenile size table:", dim(modelsuite$juv_size_table)[1]))
  } else if (!is.logical(modelsuite$juv_size_model)) {
    writeLines("\nNumber of models in juvenile size table: 1")
  } else {
    writeLines("\nJuvenile size table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_reproduction_model) & is.element("model.selection", class(modelsuite$juv_reproduction_table))) {
    writeLines(paste0("\nNumber of models in juvenile reproduction table:", dim(modelsuite$juv_reproduction_table)[1]))
  } else if (!is.logical(modelsuite$juv_reproduction_model)) {
    writeLines("\nNumber of models in juvenile reproduction table: 1")
  } else {
    writeLines("\nJuvenile reproduction table not estimated")
  }
  
  
  writeLines("\n\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nGeneral model parameter names (column 1), and specific names used in these models (column 2):")
  print.data.frame(modelsuite$paramnames[c(1,2)])
  
  writeLines("\n\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nQuality control:\n")
  if (modelsuite$qc[1,2] > 0) {
    writeLines(paste0("Survival estimated with ", modelsuite$qc[1,2], " individuals and ", modelsuite$qc[1,3], " individual transitions."))
  } else {
    writeLines("Survival not estimated.")
  }
  
  if (modelsuite$qc[2,2] > 0) {
    writeLines(paste0("Observation estimated with ", modelsuite$qc[2,2], " individuals and ", modelsuite$qc[2,3], " individual transitions."))
  } else {
    writeLines("Observation probability not estimated.")
  }
  
  if (modelsuite$qc[3,2] > 0) {
    writeLines(paste0("Size estimated with ", modelsuite$qc[3,2], " individuals and ", modelsuite$qc[3,3], " individual transitions."))
  } else {
    writeLines("Size transition not estimated.")
  }
  
  if (modelsuite$qc[4,2] > 0) {
    writeLines(paste0("Reproductive status estimated with ", modelsuite$qc[4,2], " individuals and ", modelsuite$qc[4,3], " individual transitions."))
  } else {
    writeLines("Reproduction probability not estimated.")
  }
  
  if (modelsuite$qc[5,2] > 0) {
    writeLines(paste0("Fecundity estimated with ", modelsuite$qc[5,2], " individuals and ", modelsuite$qc[5,3], " individual transitions."))
  } else {
    writeLines("Fecundity not estimated.")
  }
  
  if (modelsuite$qc[6,2] > 0) {
    writeLines(paste0("Juvenile survival estimated with ", modelsuite$qc[6,2], " individuals and ", modelsuite$qc[6,3], " individual transitions."))
  } else {
    writeLines("Juvenile survival not estimated.")
  }
  
  if (modelsuite$qc[7,2] > 0) {
    writeLines(paste0("Juvenile observation estimated with ", modelsuite$qc[7,2], " individuals and ", modelsuite$qc[7,3], " individual transitions."))
  } else {
    writeLines("Juvenile observation probability not estimated.")
  }
  
  if (modelsuite$qc[8,2] > 0) {
    writeLines(paste0("Juvenile size estimated with ", modelsuite$qc[8,2], " individuals and ", modelsuite$qc[8,3], " individual transitions."))
  } else {
    writeLines("Juvenile size transition not estimated.")
  }
  
  if (modelsuite$qc[9,2] > 0) {
    writeLines(paste0("Juvenile reproduction estimated with ", modelsuite$qc[9,2], " individuals and ", modelsuite$qc[9,3], " individual transitions."))
  } else {
    writeLines("Juvenile reproduction probability not estimated.")
  }
  
  return()
}

#' Extract Required Coefficient Values From Vital Rate Models
#' 
#' \code{.modelextract()} is an S3 generic function designed to extract the
#' estimated coefficient values from linear models created to estimate vital rates
#' in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through functions supported by \code{'lefko3'}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function typically returns a list with elements corresponding to
#' different classes of coefficients.
#' 
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract <- function(model, ...) UseMethod(".modelextract")

#' Extract Required Coefficient Values From glmmTMB-estimated Vital Rate Models
#' 
#' \code{.modelextract.glmmTMB()} extracts coefficient values from linear models 
#' estimated through the \code{\link[glmmTMB]{glmmTMB}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through \code{\link[glmmTMB]{glmmTMB}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time coefficients, typically random.}
#' \item{patches}{Vector of patch coefficients, typically random.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation or disperson parameter of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.glmmTMB <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  inda2.coef <- 0
  indb2.coef <- 0
  indc2.coef <- 0
  
  inda1.coef <- 0
  indb1.coef <- 0
  indc1.coef <- 0
  
  inda2.size2.coef <- 0
  indb2.size2.coef <- 0
  indc2.size2.coef <- 0
  inda2.flw2.coef <- 0
  indb2.flw2.coef <- 0
  indc2.flw2.coef <- 0
  
  inda1.size1.coef <- 0
  indb1.size1.coef <- 0
  indc1.size1.coef <- 0
  inda1.flw1.coef <- 0
  indb1.flw1.coef <- 0
  indc1.flw1.coef <- 0
  
  inda2.indb2.coef <- 0
  inda2.indc2.coef <- 0
  indb2.indc2.coef <- 0
  inda1.indb1.coef <- 0
  inda1.indc1.coef <- 0
  indb1.indc1.coef <- 0
  
  inda2.indb1.coef <- 0
  inda1.indb2.coef <- 0
  inda2.indc1.coef <- 0
  inda1.indc2.coef <- 0
  indb2.indc1.coef <- 0
  indb1.indc2.coef <- 0
  
  yintercept <- glmmTMB::fixef(model)[["cond"]]["(Intercept)"]
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #These are zero-inflated
  size2.zi <- 0
  size1.zi <- 0
  flw2.zi <- 0
  flw1.zi <- 0
  
  size1.size2.zi <- 0
  flw1.flw2.zi <- 0
  size1.flw1.zi <- 0
  size2.flw2.zi <- 0
  size1.flw2.zi <- 0
  size2.flw1.zi <- 0
  
  age.zi <- 0
  age.size1.zi <- 0
  age.size2.zi <- 0
  age.flw1.zi <- 0
  age.flw2.zi <- 0
  
  inda2.zi <- 0
  indb2.zi <- 0
  indc2.zi <- 0
  
  inda1.zi <- 0
  indb1.zi <- 0
  indc1.zi <- 0
  
  inda2.size2.zi <- 0
  indb2.size2.zi <- 0
  indc2.size2.zi <- 0
  inda2.flw2.zi <- 0
  indb2.flw2.zi <- 0
  indc2.flw2.zi <- 0
  
  inda1.size1.zi <- 0
  indb1.size1.zi <- 0
  indc1.size1.zi <- 0
  inda1.flw1.zi <- 0
  indb1.flw1.zi <- 0
  indc1.flw1.zi <- 0
  
  inda2.indb2.zi <- 0
  inda2.indc2.zi <- 0
  indb2.indc2.zi <- 0
  inda1.indb1.zi <- 0
  inda1.indc1.zi <- 0
  indb1.indc1.zi <- 0
  
  inda2.indb1.zi <- 0
  inda1.indb2.zi <- 0
  inda2.indc1.zi <- 0
  inda1.indc2.zi <- 0
  indb2.indc1.zi <- 0
  indb1.indc2.zi <- 0
  
  yintercept.zi <- glmmTMB::fixef(model)[["zi"]]["(Intercept)"]
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.zi <- glmmTMB::fixef(model)[["zi"]][paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  
  #Here are the time and location bits
  year.coefs <- glmmTMB::ranef(model)[["cond"]][[paramnames[(which(paramnames$mainparams == "year2")), "modelparams"]]]
  if (!all(is.na(glmmTMB::ranef(model)[["zi"]][[paramnames[(which(paramnames$mainparams == "year2")), "modelparams"]]]))) {
    year.zi <- glmmTMB::ranef(model)[["zi"]][[paramnames[(which(paramnames$mainparams == "year2")), "modelparams"]]]
  } else {
    year.zi <- 0
  }
  
  if (!all(is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]))) {
    patch.coefs <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]
  } else if (!all(is.na(glmmTMB::ranef(model)[["cond"]][[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]))) {
    patch.coefs <- glmmTMB::ranef(model)[["cond"]][[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]
  } else {
    patch.coefs <- 0
  }
  
  if (!all(is.na(glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]))) {
    patch.zi <- glmmTMB::fixef(model)[["zi"]][paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]
  } else if (!all(is.na(glmmTMB::ranef(model)[["zi"]][[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]))) {
    patch.zi <- glmmTMB::ranef(model)[["zi"]][[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]
  } else {
    patch.zi <- 0
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
                indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef, inda2.flw2.coef,
                indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef, indb1.size1.coef, indc1.size1.coef,
                inda1.flw1.coef, indb1.flw1.coef, indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, 
                indb2.indc2.coef, inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
                inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef, indb1.indc2.coef, 
                
                yintercept.zi, flw1.zi, flw2.zi, size1.zi, size2.zi, flw1.flw2.zi, size1.size2.zi, 
                size1.flw1.zi, size2.flw2.zi, size2.flw1.zi, size1.flw2.zi, age.zi, age.size1.zi, 
                age.size2.zi, age.flw1.zi, age.flw2.zi, inda2.zi, indb2.zi, indc2.zi, inda1.zi,
                indb1.zi, indc1.zi, inda2.size2.zi, indb2.size2.zi, indc2.size2.zi, inda2.flw2.zi,
                indb2.flw2.zi, indc2.flw2.zi, inda1.size1.zi, indb1.size1.zi, indc1.size1.zi,
                inda1.flw1.zi, indb1.flw1.zi, indc1.flw1.zi, inda2.indb2.zi, inda2.indc2.zi, 
                indb2.indc2.zi, inda1.indb1.zi, inda1.indc1.zi, indb1.indc1.zi, inda2.indb1.zi, 
                inda1.indb2.zi, inda2.indc1.zi, inda1.indc2.zi, indb2.indc1.zi, indb1.indc2.zi)
  
  rvars <- NA
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs, patches = patch.coefs, variances = rvars, 
                    family = stats::family(model)$family, sigma = glmmTMB::sigma(model), zeroyear = year.zi, 
                    zeropatch = patch.zi)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    if (year.as.random == TRUE) {
      newdevs <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$years[,1]))
      newdevs.zi <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$zeroyear[,1]))
    } else {
      newdevs <- rep(0, (length(yeardiff)))
      newdevs.zi <- rep(0, (length(yeardiff)))
    }
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
    
    coef.list$zeroyear <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs.zi, ncol = 1)), 
                                                    names(coef.list$zeroyear)), coef.list$zeroyear)
    yearlabels.zi <- c(yeardiff, rownames(coef.list$zeroyear)[(length(yeardiff) + 1):(length(coef.list$zeroyear[,1]))])
    coef.list$zeroyear <- as.data.frame(coef.list$zeroyear[order(yearlabels.zi),])
    rownames(coef.list$zeroyear) <- sort(yearlabels.zi)
  }
  
  coef.list$years <- coef.list$years[,1]
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      if (patch.as.random == TRUE) {
        newdevs <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$patches[,1]))
        newdevs.zi <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$zeropatch[,1]))
      } else {
        newdevs <- rep(0, (length(patchdiff)))
        newdevs.zi <- rep(0, (length(patchdiff)))
      }
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), 
                                            coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
      
      coef.list$zeropatch <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs.zi, ncol = 1)), names(coef.list$zeropatch)), 
                                              coef.list$zeropatch)
      patchlabels.zi <- c(patchdiff, rownames(coef.list$zeropatch)[(length(patchdiff) + 1):(length(coef.list$zeropatch[,1]))])
      coef.list$zeropatch <- as.data.frame(coef.list$zeropatch[order(patchlabels.zi),])
      rownames(coef.list$zeropatch) <- sort(patchlabels.zi)
    }
    coef.list$patches <- coef.list$patches[,1]
  }
  
  if (class(coef.list$zeroyear) == "data.frame") coef.list$zeroyear <- coef.list$zeroyear[,1]
  if (class(coef.list$zeropatch) == "data.frame") coef.list$zeropatch <- coef.list$zeropatch[,1]
  
  return(coef.list)
}

#' Extract Required Coefficient Values From glmer-estimated Vital Rate Models
#' 
#' \code{.modelextract.merMod()} extracts coefficient values from linear models 
#' estimated through the \code{\link[lme4]{glmer}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through \code{\link[lme4]{glmer}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time coefficients, typically random.}
#' \item{patches}{Vector of patch coefficients, typically random.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.merMod <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  inda2.coef <- 0
  indb2.coef <- 0
  indc2.coef <- 0
  
  inda1.coef <- 0
  indb1.coef <- 0
  indc1.coef <- 0
  
  inda2.size2.coef <- 0
  indb2.size2.coef <- 0
  indc2.size2.coef <- 0
  inda2.flw2.coef <- 0
  indb2.flw2.coef <- 0
  indc2.flw2.coef <- 0
  
  inda1.size1.coef <- 0
  indb1.size1.coef <- 0
  indc1.size1.coef <- 0
  inda1.flw1.coef <- 0
  indb1.flw1.coef <- 0
  indc1.flw1.coef <- 0
  
  inda2.indb2.coef <- 0
  inda2.indc2.coef <- 0
  indb2.indc2.coef <- 0
  inda1.indb1.coef <- 0
  inda1.indc1.coef <- 0
  indb1.indc1.coef <- 0
  
  inda2.indb1.coef <- 0
  inda1.indb2.coef <- 0
  inda2.indc1.coef <- 0
  inda1.indc2.coef <- 0
  indb2.indc1.coef <- 0
  indb1.indc2.coef <- 0
  
  yintercept <- lme4::fixef(model)["(Intercept)"]
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #These are the final bits related to time and location
  year.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "year2")), "modelparams"]]]
  
  if (!all(is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]))) {
    patch.coefs <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]
  } else if (!all(is.na(lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]))) {
    patch.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]
  } else {
    patch.coefs <- NA
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
                indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef, inda2.flw2.coef,
                indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef, indb1.size1.coef, indc1.size1.coef,
                inda1.flw1.coef, indb1.flw1.coef, indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, 
                indb2.indc2.coef, inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
                inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef, indb1.indc2.coef, 
                rep(0, 46))
  
  rvars <- as.data.frame(lme4::VarCorr(model))
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs, patches = patch.coefs, variances = rvars, 
                    family = stats::family(model)$family, sigma = NA, zeroyear = NA, zeropatch = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    if (year.as.random == TRUE) {
      newdevs <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$years[,1]))
    } else {
      newdevs <- rep(0, (length(yeardiff)))
    }
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      if (patch.as.random == TRUE) {
        newdevs <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$patches[,1]))
      } else {
        newdevs <- rep(0, (length(patchdiff)))
      }
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Extract Required Coefficient Values From lmer-estimated Vital Rate Models
#' 
#' \code{.modelextract.lmerMod()} extracts coefficient values from linear models 
#' estimated through the \code{\link[lme4]{lmer}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through \code{\link[lme4]{lmer}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time coefficients, typically random.}
#' \item{patches}{Vector of patch coefficients, typically random.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.lmerMod <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  inda2.coef <- 0
  indb2.coef <- 0
  indc2.coef <- 0
  
  inda1.coef <- 0
  indb1.coef <- 0
  indc1.coef <- 0
  
  inda2.size2.coef <- 0
  indb2.size2.coef <- 0
  indc2.size2.coef <- 0
  inda2.flw2.coef <- 0
  indb2.flw2.coef <- 0
  indc2.flw2.coef <- 0
  
  inda1.size1.coef <- 0
  indb1.size1.coef <- 0
  indc1.size1.coef <- 0
  inda1.flw1.coef <- 0
  indb1.flw1.coef <- 0
  indc1.flw1.coef <- 0
  
  inda2.indb2.coef <- 0
  inda2.indc2.coef <- 0
  indb2.indc2.coef <- 0
  inda1.indb1.coef <- 0
  inda1.indc1.coef <- 0
  indb1.indc1.coef <- 0
  
  inda2.indb1.coef <- 0
  inda1.indb2.coef <- 0
  inda2.indc1.coef <- 0
  inda1.indc2.coef <- 0
  indb2.indc1.coef <- 0
  indb1.indc2.coef <- 0
  
  yintercept <- lme4::fixef(model)["(Intercept)"]
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #Here are the time and location bits
  year.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "year2")), "modelparams"]]]
  
  if (!all(is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]))) {
    patch.coefs <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]
  } else if (!all(is.na(lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]))) {
    patch.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), "modelparams"]]]
  } else {
    patch.coefs <- NA
  }
  
  sigmax <- attr(summary(model)$varcor, "sc")
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
                indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef, inda2.flw2.coef,
                indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef, indb1.size1.coef, indc1.size1.coef,
                inda1.flw1.coef, indb1.flw1.coef, indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, 
                indb2.indc2.coef, inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
                inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef, indb1.indc2.coef, 
                rep(0, 46))
  
  rvars <- as.data.frame(lme4::VarCorr(model))
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs, patches = patch.coefs, variances = rvars, 
                    family = "gaussian", sigma = sigmax, zeroyear = NA, zeropatch = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    if (year.as.random == TRUE) {
      newdevs <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$years[,1]))
    } else {
      newdevs <- rep(0, (length(yeardiff)))
    }
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      if (patch.as.random == TRUE) {
        newdevs <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$patches[,1]))
      } else {
        newdevs <- rep(0, (length(patchdiff)))
      }
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Extract Required Coefficient Values From glmmTMB-estimated Vital Rate Models
#' 
#' \code{.modelextract.zeroinfl()} extracts coefficient values from linear models 
#' estimated through the \code{\link[pscl]{zeroinfl}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()}, 
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through \code{\link[glmmTMB]{glmmTMB}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time coefficients, typically random.}
#' \item{patches}{Vector of patch coefficients, typically random.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.zeroinfl <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  family.coef <- NA
  
  if (model$dist == "binomial") {family.coef <- "binomial"}
  if (model$dist == "poisson") {family.coef <- "poisson"}
  if (model$dist == "negbin") {family.coef <- "negbin"}
  
  if (is.na(family.coef)) {stop("Model distribution not recognized.")}
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  inda2.coef <- 0
  indb2.coef <- 0
  indc2.coef <- 0
  
  inda1.coef <- 0
  indb1.coef <- 0
  indc1.coef <- 0
  
  inda2.size2.coef <- 0
  indb2.size2.coef <- 0
  indc2.size2.coef <- 0
  inda2.flw2.coef <- 0
  indb2.flw2.coef <- 0
  indc2.flw2.coef <- 0
  
  inda1.size1.coef <- 0
  indb1.size1.coef <- 0
  indc1.size1.coef <- 0
  inda1.flw1.coef <- 0
  indb1.flw1.coef <- 0
  indc1.flw1.coef <- 0
  
  inda2.indb2.coef <- 0
  inda2.indc2.coef <- 0
  indb2.indc2.coef <- 0
  inda1.indb1.coef <- 0
  inda1.indc1.coef <- 0
  indb1.indc1.coef <- 0
  
  inda2.indb1.coef <- 0
  inda1.indb2.coef <- 0
  inda2.indc1.coef <- 0
  inda1.indc2.coef <- 0
  indb2.indc1.coef <- 0
  indb1.indc2.coef <- 0
  
  yintercept <- model$coefficients$count["(Intercept)"]
  
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.coef <- model$coefficients$count[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.coef <- model$coefficients$count[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #These are zero-inflated
  size2.zi <- 0
  size1.zi <- 0
  flw2.zi <- 0
  flw1.zi <- 0
  
  size1.size2.zi <- 0
  flw1.flw2.zi <- 0
  size1.flw1.zi <- 0
  size2.flw2.zi <- 0
  size1.flw2.zi <- 0
  size2.flw1.zi <- 0
  
  age.zi <- 0
  age.size1.zi <- 0
  age.size2.zi <- 0
  age.flw1.zi <- 0
  age.flw2.zi <- 0
  
  inda2.zi <- 0
  indb2.zi <- 0
  indc2.zi <- 0
  
  inda1.zi <- 0
  indb1.zi <- 0
  indc1.zi <- 0
  
  inda2.size2.zi <- 0
  indb2.size2.zi <- 0
  indc2.size2.zi <- 0
  inda2.flw2.zi <- 0
  indb2.flw2.zi <- 0
  indc2.flw2.zi <- 0
  
  inda1.size1.zi <- 0
  indb1.size1.zi <- 0
  indc1.size1.zi <- 0
  inda1.flw1.zi <- 0
  indb1.flw1.zi <- 0
  indc1.flw1.zi <- 0
  
  inda2.indb2.zi <- 0
  inda2.indc2.zi <- 0
  indb2.indc2.zi <- 0
  inda1.indb1.zi <- 0
  inda1.indc1.zi <- 0
  indb1.indc1.zi <- 0
  
  inda2.indb1.zi <- 0
  inda1.indb2.zi <- 0
  inda2.indc1.zi <- 0
  inda1.indc2.zi <- 0
  indb2.indc1.zi <- 0
  indb1.indc2.zi <- 0
  
  yintercept.zi <- model$coefficients$zero["(Intercept)"]
  
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.zi <- model$coefficients$zero[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.zi <- model$coefficients$zero[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #These are the time and location bits
  year.coef <- 0
  year.trace <- apply(as.matrix(names(model$coefficients$count)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)})
  
  if (sum(year.trace) > 0) {
    year.vec <- rep(0, (length(which(year.trace == TRUE)) + 1))
    year.vec[2:length(year.vec)] <- model$coefficients$count[which(year.trace == TRUE)]
    year.name.vec <- c(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], names(model$coefficients$count)[which(apply(as.matrix(names(model$coefficients$count)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)}) == TRUE)])
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], "", X)})
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub("as.factor()", "", X, fixed = TRUE)})
    names(year.vec) <- year.name.vec
  } else {
    if (length(mainyears) > 0) {
      year.vec <- rep(0, length(mainyears))
      names(year.vec) <- mainyears
    } else {
      year.vec <- 0
    }
    if (year.as.random) {warning("Year was not included in count modeling. Consider setting year.as.random to FALSE.")}
  }
  
  year.coef.zi <- 0
  year.trace.zi <- apply(as.matrix(names(model$coefficients$zero)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)})
  
  if (sum(year.trace.zi) > 0) {
    year.vec.zi <- rep(0, (length(which(year.trace.zi == TRUE)) + 1))
    year.vec.zi[2:length(year.vec.zi)] <- model$coefficients$zero[which(year.trace.zi == TRUE)]
    year.name.vec <- c(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], names(model$coefficients$zero)[which(apply(as.matrix(names(model$coefficients$zero)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)}) == TRUE)])
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], "", X)})
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub("as.factor()", "", X, fixed = TRUE)})
    names(year.vec.zi) <- year.name.vec
  } else {
    if (length(mainyears) > 0) {
      year.vec.zi <- rep(0, length(mainyears))
      names(year.vec.zi) <- mainyears
    } else {
      year.vec.zi <- 0
    }
    if (year.as.random) {warning("Year was not included in zero inflation modeling. Consider setting year.as.random to FALSE.")}
  }
  
  patch.vec <- NA
  patch.trace <- apply(as.matrix(names(model$coefficients$count)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)})
  
  if (sum(patch.trace) > 0) {
    patch.vec <- rep(0, (length(which(patch.trace == TRUE)) + 1))
    patch.vec[2:length(patch.vec)] <- model$coefficients$count[which(patch.trace == TRUE)]
    patch.name.vec <- c(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], names(model$coefficients$count)[which(apply(as.matrix(names(model$coefficients$count)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)}) == TRUE)])
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], "", X)})
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub("as.factor\\()", "", X)})
    names(patch.vec) <- patch.name.vec
  }
  
  patch.vec.zi <- NA
  patch.trace.zi <- apply(as.matrix(names(model$coefficients$zero)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)})
  
  if (sum(patch.trace.zi) > 0) {
    patch.vec.zi <- rep(0, (length(which(patch.trace.zi == TRUE)) + 1))
    patch.vec.zi[2:length(patch.vec.zi)] <- model$coefficients$zero[which(patch.trace.zi == TRUE)]
    patch.name.vec <- c(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], names(model$coefficients$zero)[which(apply(as.matrix(names(model$coefficients$zero)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)}) == TRUE)])
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], "", X)})
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub("as.factor\\()", "", X)})
    names(patch.vec.zi) <- patch.name.vec
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
                indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef, inda2.flw2.coef,
                indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef, indb1.size1.coef, indc1.size1.coef,
                inda1.flw1.coef, indb1.flw1.coef, indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, 
                indb2.indc2.coef, inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
                inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef, indb1.indc2.coef, 
                
                yintercept.zi, flw1.zi, flw2.zi, size1.zi, size2.zi, flw1.flw2.zi, size1.size2.zi, 
                size1.flw1.zi, size2.flw2.zi, size2.flw1.zi, size1.flw2.zi, age.zi, age.size1.zi, 
                age.size2.zi, age.flw1.zi, age.flw2.zi, inda2.zi, indb2.zi, indc2.zi, inda1.zi,
                indb1.zi, indc1.zi, inda2.size2.zi, indb2.size2.zi, indc2.size2.zi, inda2.flw2.zi,
                indb2.flw2.zi, indc2.flw2.zi, inda1.size1.zi, indb1.size1.zi, indc1.size1.zi,
                inda1.flw1.zi, indb1.flw1.zi, indc1.flw1.zi, inda2.indb2.zi, inda2.indc2.zi, 
                indb2.indc2.zi, inda1.indb1.zi, inda1.indc1.zi, indb1.indc1.zi, inda2.indb1.zi, 
                inda1.indb2.zi, inda2.indc1.zi, inda1.indc2.zi, indb2.indc1.zi, indb1.indc2.zi)
  
  rvars <- NA
  
  if (family.coef == "negbin") {
    sigma <- model$theta
  } else {
    sigma <- 1
  }
  
  coef.list <- list(coefficients = coef.vec, years = year.vec, patches = patch.vec, variances = rvars, 
                    family = model$dist, sigma = sigma, zeroyear = year.vec.zi, zeropatch = patch.vec.zi)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, names(coef.list$years))
  
  if (length(yeardiff) > 0) {
    
    names(coef.list$years)[which(names(coef.list$years) == "")] <- yeardiff
    
  }
  
  yeardiff <- setdiff(mainyears, names(coef.list$zeroyear))
  
  if (length(yeardiff) > 0) {
    
    names(coef.list$zeroyear)[which(names(coef.list$zeroyear) == "")] <- yeardiff
    
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      names(coef.list$patches)[which(names(coef.list$patches) == "")] <- patchdiff
      
    }
  }
  
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$zeropatch))
    
    if (length(patchdiff) > 0) {
      names(coef.list$zeropatch)[which(names(coef.list$zeropatch) == "")] <- patchdiff
      
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  if (class(coef.list$zeroyear) == "data.frame") coef.list$zeroyear <- coef.list$zeroyear[,1]
  if (class(coef.list$zeropatch) == "data.frame") coef.list$zeropatch <- coef.list$zeropatch[,1]
  
  return(coef.list)
}

#' Extract Required Coefficient Values From glm-estimated Vital Rate Models
#' 
#' \code{.modelextract.glm()} extracts coefficient values from linear models 
#' estimated through the \code{\link[stats]{glm}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through \code{\link[stats]{glm}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time coefficients.}
#' \item{patches}{Vector of patch coefficients.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.glm <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  family.coef <- NA
  
  if (model$family["family"] == "binomial") {family.coef <- "binomial"}
  if (model$family["family"] == "poisson") {family.coef <- "poisson"}
  if (is.element("negbin", class(model))) {family.coef <- "negbin"}
  
  if (is.na(family.coef)) {stop("Model distribution not recognized.")}
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  inda2.coef <- 0
  indb2.coef <- 0
  indc2.coef <- 0
  
  inda1.coef <- 0
  indb1.coef <- 0
  indc1.coef <- 0
  
  inda2.size2.coef <- 0
  indb2.size2.coef <- 0
  indc2.size2.coef <- 0
  inda2.flw2.coef <- 0
  indb2.flw2.coef <- 0
  indc2.flw2.coef <- 0
  
  inda1.size1.coef <- 0
  indb1.size1.coef <- 0
  indc1.size1.coef <- 0
  inda1.flw1.coef <- 0
  indb1.flw1.coef <- 0
  indc1.flw1.coef <- 0
  
  inda2.indb2.coef <- 0
  inda2.indc2.coef <- 0
  indb2.indc2.coef <- 0
  inda1.indb1.coef <- 0
  inda1.indc1.coef <- 0
  indb1.indc1.coef <- 0
  
  inda2.indb1.coef <- 0
  inda1.indb2.coef <- 0
  inda2.indc1.coef <- 0
  inda1.indc2.coef <- 0
  indb2.indc1.coef <- 0
  indb1.indc2.coef <- 0
  
  yintercept <- model$coefficients["(Intercept)"]
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #These are the time and location bits
  year.coef <- 0
  year.trace <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)})
  
  if (sum(year.trace) > 0) {
    year.vec <- rep(0, (length(which(year.trace == TRUE)) + 1))
    year.vec[2:length(year.vec)] <- model$coefficients[which(year.trace == TRUE)]
    year.name.vec <- c(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)}) == TRUE)])
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], "", X)})
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub("as.factor()", "", X, fixed = TRUE)})
    names(year.vec) <- year.name.vec
  } else {
    if (length(mainyears) > 0) {
      year.vec <- rep(0, length(mainyears))
      names(year.vec) <- mainyears
    } else {
      year.vec <- 0
    }
    if (year.as.random) {warning("Year was not included in the modeling. Consider setting year.as.random to FALSE.")}
  }
  
  patch.vec <- NA
  patch.trace <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)})
  
  if (sum(patch.trace) > 0) {
    patch.vec <- rep(0, (length(which(patch.trace == TRUE)) + 1))
    patch.vec[2:length(patch.vec)] <- model$coefficients[which(patch.trace == TRUE)]
    patch.name.vec <- c(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)}) == TRUE)])
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], "", X)})
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub("as.factor\\()", "", X)})
    names(patch.vec) <- patch.name.vec
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
                indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef, inda2.flw2.coef,
                indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef, indb1.size1.coef, indc1.size1.coef,
                inda1.flw1.coef, indb1.flw1.coef, indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, 
                indb2.indc2.coef, inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
                inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef, indb1.indc2.coef, 
                rep(0, 46))
  
  if (is.element("negbin", class(model))) {
    sigma <- model$theta
  } else if (is.element("gaussian", class(model))) {
    sigma <- stats::sigma(model)
  } else {
    sigma <- NA
  }
  
  coef.list <- list(coefficients = coef.vec, years = year.vec, patches = patch.vec, variances = NA, 
                    family = family.coef, sigma = sigma, zeroyear = NA, zeropatch = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, names(coef.list$years))
  
  if (length(yeardiff) > 0) {
    
    names(coef.list$years)[which(names(coef.list$years) == "")] <- yeardiff
    
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      names(coef.list$patches)[which(names(coef.list$patches) == "")] <- patchdiff
      
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Extract Required Coefficient Values From lm-estimated Vital Rate Models
#' 
#' \code{.modelextract.lm()} extracts coefficient values from linear models 
#' estimated through the \code{\link[stats]{lm}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through \code{\link[stats]{lm}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time coefficients.}
#' \item{patches}{Vector of patch coefficients.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.lm <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  inda2.coef <- 0
  indb2.coef <- 0
  indc2.coef <- 0
  
  inda1.coef <- 0
  indb1.coef <- 0
  indc1.coef <- 0
  
  inda2.size2.coef <- 0
  indb2.size2.coef <- 0
  indc2.size2.coef <- 0
  inda2.flw2.coef <- 0
  indb2.flw2.coef <- 0
  indc2.flw2.coef <- 0
  
  inda1.size1.coef <- 0
  indb1.size1.coef <- 0
  indc1.size1.coef <- 0
  inda1.flw1.coef <- 0
  indb1.flw1.coef <- 0
  indc1.flw1.coef <- 0
  
  inda2.indb2.coef <- 0
  inda2.indc2.coef <- 0
  indb2.indc2.coef <- 0
  inda1.indb1.coef <- 0
  inda1.indc1.coef <- 0
  indb1.indc1.coef <- 0
  
  inda2.indb1.coef <- 0
  inda1.indb2.coef <- 0
  inda2.indc1.coef <- 0
  inda1.indc2.coef <- 0
  indb2.indc1.coef <- 0
  indb1.indc2.coef <- 0
  
  yintercept <- model$coefficients["(Intercept)"]
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]])) {flw2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]])) {flw1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]])) {size2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]])) {size1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]])) {age.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "age")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "age")), "modelparams"])]}
  
  #Here are the individual covariate bits
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]])) {inda2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]])) {indb2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]])) {indc2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]])) {inda1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]])) {indb1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]])) {indc1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {inda2.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indb2.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])])) {indc2.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {inda2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indb2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])])) {indc2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {inda1.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indb1.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])])) {indc1.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "size1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {inda1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indb1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])])) {indc1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "repst1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda2.indb2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda2.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb2.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda1.indb1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda1.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb1.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])])) {inda2.indb1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])])) {inda1.indb2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {inda2.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {inda1.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcova1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])])) {indb2.indc1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb2")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc1")), "modelparams"])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])])) {indb1.indc2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "indcovb1")), "modelparams"], ":", paramnames[(which(paramnames$mainparams == "indcovc2")), "modelparams"])]}
  
  #These are the time and location bits
  year.coef <- 0
  year.trace <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)})
  
  if (sum(year.trace) > 0) {
    year.vec <- rep(0, (length(which(year.trace == TRUE)) + 1))
    year.vec[2:length(year.vec)] <- model$coefficients[which(year.trace == TRUE)]
    year.name.vec <- c(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], X)}) == TRUE)])
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "year2")), "modelparams"], "", X)})
    names(year.vec) <- year.name.vec
  } else {
    if (length(mainyears) > 0) {
      year.vec <- rep(0, length(mainyears))
      names(year.vec) <- mainyears
    } else {
      year.vec <- 0
    }
    if (year.as.random) {warning("Year was not included in the modeling. Consider setting year.as.random to FALSE.")}
  }
  
  patch.vec <- NA
  patch.trace <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)})
  
  if (sum(patch.trace) > 0) {
    patch.vec <- rep(0, (length(which(patch.trace == TRUE)) + 1))
    patch.vec[2:length(patch.vec)] <- model$coefficients[which(patch.trace == TRUE)]
    patch.name.vec <- c(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], X)}) == TRUE)])
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "patch")), "modelparams"], "", X)})
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub("as.factor\\()", "", X)})
    names(patch.vec) <- patch.name.vec
  }
  
  sigmax <- stats::sigma(model)
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
                indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef, inda2.flw2.coef,
                indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef, indb1.size1.coef, indc1.size1.coef,
                inda1.flw1.coef, indb1.flw1.coef, indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, 
                indb2.indc2.coef, inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
                inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef, indb1.indc2.coef, 
                rep(0, 46))
  
  coef.list <- list(coefficients = coef.vec, years = year.vec, patches = patch.vec, variances = NA, 
                    family = "gaussian", sigma = sigmax, zeroyear = NA, zeropatch = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, names(coef.list$years))
  
  if (length(yeardiff) > 0) {
    
    names(coef.list$years)[which(names(coef.list$years) == "")] <- yeardiff
    
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      names(coef.list$patches)[which(names(coef.list$patches) == "")] <- patchdiff
      
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Return Simple List When Model Is Scalar Numeric
#' 
#' \code{.modelextract.numeric()} returns NA when a vital rate model is simply a
#' scalar numeric value. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through functions supported by \code{'lefko3'}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with single values for matrix estimation.
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.numeric <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  coef.list <- list(coefficients = c(model), years = c(0), patches = c(0), variances = as.data.frame(0), 
                    family = 1, sigma = 1, zeroyear = NA, zeropatch = NA)
  
  return(coef.list)
}

#' Return Simple List When Model Is Logical
#' 
#' \code{.modelextract.logical()} returns NA when a vital rate model is simply a
#' logical value. Used to supply coefficients to \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, and \code{\link{aflefko2}()}.
#' 
#' @param model Model estimated through functions supported by \code{'lefko3'}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns NA.
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.logical <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  coef.list <- list(coefficients = c(0), years = c(0), patches = c(0), variances = as.data.frame(0), 
                    family = 1, sigma = 1, zeroyear = NA, zeropatch = NA)
  
  return(coef.list)
}
