#' @title Generate dataset of normal distributed repeated observations in a one subgroup design
#' @description \code{r.gee.1subgroup} generates data for a design with one subgroup within a full population. Each baseline-observation is normal distributed with mean \deqn{\beta_0} in placebo group and \deqn{\beta_0+\beta_1} in treatment group.
#' Measurements after baseline have mean \deqn{\beta_0+\beta_2*t} in placebo group and \deqn{\beta_0+\beta_1+\beta_2*t+\beta_3*t} in treatment group where \deqn{t} is the measurement time. Whether the effect can be found solely in the subgroup or additionally a certain amount outside of the subgroup can be specified as well as a potential different covariance-structure within subgroup and in the complementary subgroup.
#'
#' @param n        overall sample size for the overall population
#' @param reg      list containing coefficients \deqn{\beta_0} to \deqn{\beta_0} for complementary population, \code{reg[[1]]} and subpopulation, \code{reg[[2]]}: see 'Details'.
#' @param sigma    vector with standard deviations for generated observations c(complementary population, subpopulation).
#' @param rho      variable used together with \code{theta} to describe correlation between two adjacent timepoints: see 'Details'.
#' @param theta    variable used together with \code{rho} to describe correlation between two adjacent timepoints: see 'Details'.
#' @param tau      subgroup prevalence.
#' @param k        sample size allocation factor between treatment groups: see 'Details'.
#' @param Time     list of timepoints \eqn{t} that have to be generated: see 'Details'.
#' @param OD       percentage of observed overall dropout at last timepoint: see 'Details'.
#'
#' @details
#' For \code{reg}\code{list}(c(\eqn{\beta_0^F\S,\beta_1^F\S,\beta_2^F\S,\beta_3^F\S}), c(\eqn{\beta_0^S,\beta_1^S,\beta_2^S,\beta_3^S})) and variances \code{sigma}=(\eqn{\sigma_F\S, \sigma_S}) function \code{r.gee.1subgroup} generates data given correlation-variables \eqn{\rho} and \eqn{\theta} as follows (and let t=0 be the baseline measurement):
#'
#' Placebo group - complementary population \eqn{y_{it}=N(\beta_0+\beta_2*t,\sigma_F\S)},
#' Placebo group - within subgroup \eqn{y_{it}=N(\beta_0+\beta_2*t,\sigma_S)},
#' Treatment group - complementary population \eqn{y_{it}=N(\beta_0+\beta_1+\beta_2*t+\beta_3*t,\sigma_F\S)},
#' Treatment group - within subgroup \eqn{y_{it}=N(\beta_0+\beta_1+\beta_2*t+\beta_3*t,\sigma_S)}.
#' Correlation between measurements - \eqn{corr(\epsilon_it,\epsilon_io)=\rho^{(t-o)^\theta}}
#'
#'  Argument \code{k} is the sample size allocation factor, i.e. the ratio between control and treatment. Let \eqn{n_C} and \eqn{n_T} denote sample sizes of control and treatment groups respectively, then \eqn{k = n_T/n_C}.
#'
#' Argument \code{Time} is the vector denoting all measuring-times, i. e. every value for \eqn{t}.
#'
#' Argument \code{OD} sets the overall dropout rate observed at the last timepoint. For \code{OD}=0.5, 50 percent of all observation had a dropout event at some point. If a subject experienced a dropout the starting time of the dropout is equally distributed over all timepoints.
#'
#' @return \code{r.gee.1subgroup} returns a list with 7 different entries. Every Matrix rows are the simulated subjects and the columns are the observed time points.
#'
#' The first list element is a vector containing subject ids.
#' The second element contains a matrix with the outcomes of a subject with row being the subjects and columns being the measuring-timepoints
#' Elements 3 to 5 return matrices with the information of which patients have baseline-measurements, which patients belong to treatment and which to control and what are the observed timepoints for each patient respectively.
#' The sixth entry returns a matrix which contains the residuals of each measurement.
#' The seventh entry returns the sub-population identification.
#'
#' @source \code{r.gee.1subgroup} uses code contributed by Roland Gerard Gera
#'
#' @examples
#'
#' set.seed(2015)
#' dataset<-r.gee.1subgroup(n=200, reg=list(c(0,0,0,0.1),c(0,0,0,0.1)), sigma=c(3,2.5),
#' tau=0.5, rho=0.25, theta=1, k=1.5, Time=c(0:5), OD=0)
#' dataset
#' @import MASS
#' @export

r.gee.1subgroup<-function(n, reg, sigma, rho ,theta , tau, k, Time, OD){

  size_S  = round(n*tau)
  size_SC = n-size_S

  #  if no data from the subpopulation is generated
  if (size_S==0) {

    Datenliste_S = c()

  } else {

    Datenliste_S <-Datagen(Varianz = sigma[2] ,
                           rho = rho ,
                           theta = theta ,
                           k = k,
                           Koeffizienten = reg[[2]] ,
                           n = size_S ,
                           Time = Time,
                           OverallDropout = OD,
                           Typ="S")

  }
  # Falls keine Daten aus der komplement?ren Subgruppe generiert werden m?ssen
  if ( size_SC == 0 ) {

    Datenliste_SC <-c()

  } else {

    Datenliste_SC <- Datagen(Varianz = sigma[1] ,
                             rho = rho ,
                             theta = theta ,
                             Koeffizienten = reg[[1]] ,
                             k = k,
                             n = size_SC ,
                             Time = Time ,
                             OverallDropout = OD,
                             Typ="SC")
    #  Falls Daten aus S erzeugt wurden, incrementiere die ID aus SC, um diese Anzupassen
    if (size_S!=0) Datenliste_SC$id<-Datenliste_SC$id+max(Datenliste_S$id)
  }


  Datenliste <- list()
  # Fuege die Daten zusammen
  for (i in 1:length(Datenliste_SC)){

    Datenliste[[i]] <-rbind(Datenliste_S[[i]] , Datenliste_SC[[i]])

  }

  # Gebe den Listeneintr?gen die Namen aus SC
  names(Datenliste)=names(Datenliste_S)

  return(Datenliste)
}


# Function which generates the List


# Generierung von Simulationsdaten, H1 True und H1 False ------------------
Datagen <- function(Varianz, rho , theta , Koeffizienten , k , n , Time , OverallDropout , Typ){
  #require(MASS)
  kp=1/(k+1)
  TimePoints=length(Time)
  # Erstelle eine "Maske", fuer NA um dropouts zu erzeugen
  Remain.matrix=matrix(1 , nrow = n , ncol = TimePoints)

  # only if droppout is != 0
  if(OverallDropout != 0){
    # calculate which patients expereince a dropout

    drops=sample(1:n , OverallDropout*n)

    # the range in which dropouts can occure
    range=2:length(Time)

    for(i in 1:length(drops)){

      # draw where at random the missing begins
      k=sample(range,1)
      #calculate matrix where 1 is the part where no dropout occured and NA where dropout exists
      Remain.matrix[drops[i],k:length(Time)]=NA

    }
  }

  ## Bestimme die Korrelation, die zwischen den Daten bestehen soll -----
  Korrelation <- gen_cov_cor(var=Varianz,
                             rho=rho,
                             theta=theta,
                             Time=Time,
                             cov=TRUE)


  # Berechne die einzelnen Parameter der Formel -----------------------------
  id = matrix(1 , ncol=TimePoints , nrow=n) * 1:n

  error = mvrnorm(n , mu=matrix(0,ncol=TimePoints) , Sigma=Korrelation) * Remain.matrix

  Baseline = matrix(1 , nrow=n , ncol = TimePoints)

  Gr = rbind(matrix(0 , nrow = floor(n*kp)   , ncol = TimePoints),
             matrix(1 , nrow = n-floor(n*kp) , ncol = TimePoints))

  time = t( matrix(rep(Time,n) , nrow = TimePoints) )

  Gr_x_Time = time*Gr

  y = Baseline*Koeffizienten[1]+
    Gr*Koeffizienten[2]+
    time*Koeffizienten[3]+
    Gr_x_Time*Koeffizienten[4]+
    error

  population = matrix(Typ , ncol = TimePoints , nrow = n)

  # Fuege alle einzelne "Messdaten zusammen" ------------
  study = list(id = id ,
               y=y ,
               Baseline=Baseline,
               Gr=Gr,
               Time=time,
               error=error,
               population = population)
  class(study)="repdata"
  return(study)
}
