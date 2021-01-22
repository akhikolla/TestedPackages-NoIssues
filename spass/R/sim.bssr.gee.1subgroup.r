#' @title Simulation of a longitudinal one subgroup design with internal pilot Study
#'
#' @description Given estimates of the treatment effects to be proven, the variances, and the prevalence,
#' \code{sim.bssr.gee.1subgroup} calculates an initial sample size and performs a blinded sample size recalculation
#' after a pre-specified number of subjects have been enrolled. Each observation is simulated and a final analysis executed.
#' Several variations are included, such as different approximations or sample size allocation.
#'
#' @param nsim         number of simulation runs.
#' @param alpha        level (type I error) to which the hypothesis is tested.
#' @param tail         which type of test is used, e.g. which quartile und H0 is calculated
#' @param beta         type II error (power=1-beta) to which an alternative should be proven.
#' @param delta        vector of true treatment effects, c(overall population, inside subgroup).
#' @param sigma_pop    vector of true standard deviations of the treatment effects, c(overall population, subgroup).
#' @param tau          subgroup prevalence.
#' @param rho          true correlation coefficient between two adjacent timepoints
#' @param vrho         initial expectation of the correlation coefficient between two adjacent timepoints
#' @param theta        true correlation absorption coefficient if timepoints are farther apart
#' @param vtheta       expected correlation absorption coefficient if timepoints are farther apart
#' @param vdelta       vector of treatment effects to be proven, c(overall population, inside subgroup).
#' @param vsigma_pop   vector of assumed standard deviations, c(overall population, inside subgroup).
#' @param Time         vector of measured timepoints
#' @param V            working covariance matrix.
#' @param OD           overall dropout measured at last timepoint
#' @param vdropout     vector of expected dropouts per timepoint if missingness is to be expected
#' @param missingtype  true missingtype underlying the missingness
#' @param vmissingtype initial assumptions about the missingtype underlying the missingness
#' @param rec.at       blinded sample size review is performed after \code{rec.at}*\eqn{100\%} subjects of the initial sample size calculation.
#' @param k            sample size allocation factor between groups: see 'Details'.
#' @param model        which of the two often revered statistical models should be used?: see 'Details'.
#' @param seed         set seed value for the simulations to compare results.
#'
#' @details
#' This function combines sample size estimation, blinded sample size re-estimation and analysis in a design with a subgroup within a full population where we want to test for treatment effects between a control and a treatment group.
#' The required sample size for the control and treatment group to prove an existing
#' alternative \code{delta} with a specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha} is calculated prior to the study and then recalculated in an internal pilot study.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#'
#' @return \code{sim.bssr.1subgroup} returns a data.frame containing the mean and variance of recalculated sample sizes within the control group and treatment group respectively and the achieved simulated power along with all relevant parameters.
#'
#' @source \code{sim.bssr.gee.1subgroup} uses code contributed by Roland Gerard Gera.
#'
#' @seealso \code{sim.bssr.gee.1subgroup} makes use of \code{\link{n.gee.1subgroup}}, \code{\link{bssr.gee.1subgroup}}, and \code{\link{r.gee.1subgroup}}.
#'
#' @examples
#' sim.bssr.gee.1subgroup(nsim = 5,missingtype = "intermittened")
#'
#' @import geepack
#' @import mvtnorm
#' @import multcomp
#' @import utils
#' @export

sim.bssr.gee.1subgroup<-function(nsim=1000,
                                 alpha=0.05,
                                 tail="both",
                                 beta=0.2,
                                 delta=c(0.1,0.1), vdelta=c(0.1,0.1),
                                 sigma_pop=c(3,3),vsigma_pop=c(3,3),
                                 tau=0.5,
                                 rho=0.25, vrho=0.25,
                                 theta=1, vtheta=1,
                                 Time=0:5,
                                 rec.at=0.5,
                                 k=1, model=1, V=diag(rep(1,length(Time))), OD=0,
                                 vdropout=rep(0,length(Time)),
                                 missingtype="none", vmissingtype="none", seed=2015){


  set.seed(seed)
  # First we need to know our assumed variance in our regression coefficients. For that we include every information we assume that our model has.

  # set progresbar
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  # calculate to get draftet into treatment group (1 in our model)
  prob=1-(1/(1+k))
  if(model==1){
    D=c()
    for (i in 1:length(Time)) D=rbind(D,c(1,prob,Time[i],prob*Time[i]))
    correction=rbind(c(1,1,1,1),c(1,1/prob,1,1/prob),c(1,1,1,1),c(1,1/prob,1,1/prob))
    regressor=4

  }
  if(model==2){
    D=c()
    for (i in 1:length(Time)) D=rbind(D,c(1,Time[i],prob*Time[i]))
    correction=rbind(c(1,1,1),c(1,1,1),c(1,1,1/prob))
    regressor=3
  }
  reg=list(c(0,0,0,delta[1]),c(0,0,0,delta[2]))

  # calculate the initial assumptions for the Regressionvarianceparameters.
  sigma_sq_F = sandwich(yCov = gen_cov_cor(var = vsigma_pop[1],rho = vrho,theta = vtheta,Time = Time,cov = TRUE),D = D,V = V,correctionmatrix = correction,missing = vdropout,missingtype = vmissingtype)

  sigma_sq_S = sandwich(yCov = gen_cov_cor(var = vsigma_pop[2],rho = vrho,theta = vtheta,Time = Time,cov = TRUE),D = D,V = V,correctionmatrix = correction,missing = vdropout,missingtype = vmissingtype)


  # the sd of regression estimation is
  sigma_F =sqrt(sigma_sq_F[regressor,regressor])
  sigma_S =sqrt(sigma_sq_S[regressor,regressor])


  # Calculate the initial sample size
  n_initial = n.gee.1subgroup(alpha = alpha,beta = beta,delta = vdelta,sigma =c(sigma_F,sigma_S),tau = tau,k = k)$n

  # calculate samplesize wherer a BSSR should happen
  n_bssr=round(sum(n_initial)*rec.at)



  rej=0
  n_reest=c()
  # Here the loops begin
  for (repeats in 1:nsim){
    intdata_F = r.gee.1subgroup(n=n_bssr,reg = reg,sigma = sigma_pop,rho = rho,theta = theta,tau = tau,k = k,Time = Time,OD = OD)

    intdata_S = intdata_F
    for (i in 1:length(intdata_F)){
      intdata_S[[i]] = intdata_S[[i]][(intdata_S$population=="S")[,1],]
    }

    estim_F=estimcov(data = intdata_F,Time = Time)
    estim_S=estimcov(data = intdata_S,Time = Time)

    sigma_sq_F = sandwich(yCov = gen_cov_cor(var = estim_F[[1]][1],
                                             rho = estim_F[[1]][2],
                                             theta = estim_F[[1]][3],
                                             Time = Time,cov = TRUE),
                          D = D,V = V,correctionmatrix = correction,
                          missing = estim_F[[2]],missingtype = missingtype)

    sigma_sq_S = sandwich(yCov = gen_cov_cor(var = estim_S[[1]][1],
                                             rho = estim_S[[1]][2],
                                             theta = estim_S[[1]][3],
                                             Time = Time,cov = TRUE),
                          D = D,V = V,correctionmatrix = correction,
                          missing = estim_S[[2]],missingtype = missingtype)

    sigma_F_reest = sqrt(sigma_sq_F[regressor,regressor])
    sigma_S_reest = sqrt(sigma_sq_S[regressor,regressor])

    n_reest = cbind(n_reest, sum(bssr.gee.1subgroup(alpha = alpha,beta = beta,delta = vdelta,estsigma = c(sigma_F_reest,sigma_S_reest),tau = tau,k = k)$n))

    if(n_bssr<n_reest[repeats]){
      # generate new data
      restdata=r.gee.1subgroup(n=n_reest[repeats]-n_bssr,reg = reg,sigma = sigma_pop,rho = rho,theta = theta,tau = tau,k = k,Time = Time,OD = OD)

      # increase patient id, to match the already generated data
      restdata$id=restdata$id+n_bssr


      enddata=intdata_F
      # bind both datasets
      for (i in 1:length(intdata_F)){
        enddata[[i]]=rbind(intdata_F[[i]],restdata[[i]])
      }
    }

    enddata_S = enddata
    for (i in 1:length(enddata)){
      enddata_S[[i]] = enddata_S[[i]][(enddata_S$population=="S")[,1],]
    }


    geesol = simGeeHelpfunction(In.Mat.1=enddata,model=model)

    t_F=sqrt(coefficients(summary(geesol))$Wald)[regressor]
    if (coefficients(summary(geesol))$Estimate[regressor]<0) t_F=-sqrt(coefficients(summary(geesol))$Wald)[regressor]

    geesol = simGeeHelpfunction(In.Mat.1=enddata_S,model=model)

    t_S = sqrt(coefficients(summary(geesol))$Wald)[regressor]
    if (coefficients(summary(geesol))$Estimate[regressor]<0) t_S = -sqrt(coefficients(summary(geesol))$Wald)[regressor]

    Z=c(t_F,t_S)

    intertest = rbind(c(1,sqrt(tau)),c(sqrt(tau),1))
    Critval=qmvnorm(1-alpha, delta=c(0,0), corr=intertest, tail="both")$quantile

    if (Z[1]>=Critval||Z[2]>=Critval) rej=rej+1

    setTxtProgressBar(pb, repeats)
  }

  rej=rej/nsim
  ninitial=sum(n_initial)
  meanreest=mean(n_reest)
  sdreest=sd(n_reest)

  return(c( power=rej,ninitial=ninitial,meanreest=meanreest,sdreest=sdreest))
}









simGeeHelpfunction <- function(In.Mat.1,model){

  sizes = dim(In.Mat.1$y)
  if(model==1) func=y~Gr*Time
  if(model==2) func=y~Time+GrxTime

  # Calculate weigths if gee is possible ----------------------------------------
  gewicht = c()

  for (cols in 1:sizes[2]){

    # finde Herraus wie Personen per Zeitschritt rausgeworfen werden
    gewicht[cols] = 1-sum(is.na(In.Mat.1[[2]][,cols]))/length(In.Mat.1[[2]][,cols])
  }
  # Nun Berechne die Tatsaechlich verwendeten Gewichte fuer alle Zeiten
  gewicht = 1/gewicht
  # Wende alle Gewichte auf alle Personen an
  gewicht = matrix(rep(gewicht,sizes[1]),nrow=1)
  id<-NULL
  # Convert your data in the needed format ---------------------------------------------------
  study <- data.frame(id=matrix(t(In.Mat.1[[1]]),ncol=1),
                      y=matrix(t(In.Mat.1[[2]]),ncol=1),
                      Gr=matrix(t(In.Mat.1[[4]]),ncol=1),
                      Time=matrix(t(In.Mat.1[[5]]),ncol=1),
                      GrxTime=matrix(t(In.Mat.1[[5]]*In.Mat.1[[4]]),ncol=1),
                      error=matrix(t(In.Mat.1[[7]]),ncol=1),
                      gewicht = matrix(gewicht,ncol=1)
  )

  ergebnis_gee_beta<-geeglm(formula=func,
                            id=id,
                            weights=gewicht,
                            data=study,
                            family=gaussian,
                            corstr="independence"
  )
  return(ergebnis_gee_beta)
}
