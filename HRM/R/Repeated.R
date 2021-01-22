####################################################################################################################################
### Filename:    Repeated.R
### Description: Functions that process the user input; these functions are called by generic S3 methods in the file S3methods.R
###              EEG data example
###
###
####################################################################################################################################

#' EEG data of 160 subjects
#'
#' A dataset containing EEG data (Staffen et al., 2014) of 160 subjects, 4 variables are measured at ten different locations.
#'
#' The columns are as follows:
#' \itemize{
#'   \item group. Diagnostic group of the subject: Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI), Subject Cognitive Complaints (SCC+, SCC-).
#'   \item value. Measured data of a subject at a specific variable and region.
#'   \item sex. Sex of the subject: Male (M) or Female (W).
#'   \item subject. A unique identification of a subject.
#'   \item variable. The variales measured are activity, complexity, mobility and brain rate coded from 1 to 4.
#'   \item region. Frontal left/right, central left/right, temporal left/right, occipital left/right, parietal left/right coded as 1 to 10.
#'   \item dimension. Mixing variable and region together, levels range from 1 to 40.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name EEG
#' @usage data(EEG)
#' @format A data frame with 6400 rows and 7 variables.
"EEG"


#' Test for no main treatment effect, no main time effect, no simple treatment effect and no interaction between treatment and time
#'
#' @param data A list containing the data matrices of all groups. The rows are the independent subjects, these observations are assumed to be multivariate normally distributed. The columsn of all matrices need to be in the same order.
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.matrices <- function(data, alpha=0.05){

  if(!is.list(data)){
    stop("data needs to be a list containing the data matrices of all groups")
  }

  a <- length(data)
  n <- rep(0,a)
  d <- rep(0,a)
  for(i in 1:a){
    tmp <- data[[i]]
    if(!is.matrix(tmp)){
      stop("The elements of data need to be matrices.")
    }
    n[i] <- dim(tmp)[1]
    d[i] <- dim(tmp)[2]
  }

  if(mean(d) != d[1]){
    stop("The number of measurements for each group need to be the same.")
  }
  d <- d[1]

  if(d < 2 & a < 2) {
    stop("At least two measurements per subject or two groups are needed.")
  }


  # choose appropriate tests based on the number of factors present
  testing <- rep(0, 5)

  if(a == 1 & d > 1) { testing <- c(0,0,1,0,0) }
  if(a > 1 & d == 1) { testing <- c(1,1,0,0,0) }
  if(a > 1 & d > 1) { testing <- c(1,1,1,1,1) }

  temp0 <- if(testing[1]) { hrm.A.weighted(n,a,d,data,alpha) }
  temp1 <- if(testing[2]) { hrm.A.unweighted(n,a,d,data,alpha) }
  temp2 <- if(testing[3]) { hrm.B(n,a,d,data,alpha) }
  temp3 <- if(testing[4]) { hrm.AB(n,a,d,data,alpha) }
  temp4 <- if(testing[5]) { hrm.A_B(n,a,d,data,alpha) }

  output <- list()
  output$result <- rbind(temp0,temp1,temp2,temp3,temp4)
  output$formula <- NULL
  output$alpha <- alpha
  output$subject <- NULL
  output$factors <- list(NULL, NULL)
  output$data <- data
  output$nonparametric <- FALSE
  class(output) <- "HRM"

  return (output)
}

#' Test for no main effects and interactino effects of one between-subject factor and one crossed within-subject factors
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.one <- function(X, alpha, group , factor1, subject, data, testing = rep(1,4), formula, nonparametric, np.correction ){

  ranked <- NULL
  varQGlobal <- NULL
  means <- NULL
  correction <- NULL
  tmpQ1g <- NULL
  tmpQ2g <- NULL

  temp0 <- if(testing[1]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "A", paste(as.character(group), " (weighted)"), nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g)}
  temp1 <- if(testing[2]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "Au", paste(as.character(group)), nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g)}
  temp2 <- if(testing[3]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "B", paste(as.character(factor1)), nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g)}
  temp3 <- if(testing[4]) {hrm.1w.1f(X, alpha, group , factor1, subject, data, "AB", paste(as.character(group), ":",as.character(factor1)), nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g)}
  temp4 <- if(testing[3] & (testing[1] | testing[2])) {hrm.1w.1f(X, alpha, group , factor1, subject, data, "A|B", paste(as.character(group), "|",as.character(factor1)), nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g)}
  #temp5 <- if(testing[3] & (testing[1] | testing[2])) {hrm.1w.1f(X, alpha, group , factor1, subject, data, "B|A", paste(as.character(factor1), "|",as.character(group)), nonparametric, ranked, varQGlobal, np.correction)}

  output <- list()
  output$result <- rbind(temp0, temp1, temp2, temp3, temp4)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c(factor1))
  output$data <- X
  output$mean <- means
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- correction
  class(output) <- "HRM"
  rownames(output$result) <- 1:dim(output$result)[1]
  return (output)
}

#' Test for no main effects and interactino effects of one between-subject factor and two crossed within-subject factors
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.within <- function(X, alpha, group , factor1, factor2, subject, data, testing = rep(1,7), formula, nonparametric, np.correction ){

  ranked <- NULL
  varQGlobal <- NULL
  correction <- NULL

  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- vector("list", length = 7)
  for(i in 1:7){
    if(testing[i]) {
      temp[[i]] <- hrm.1w.2f(X, alpha, group , factor1, factor2, subject, data, H = i, "", nonparametric, ranked, varQGlobal, np.correction )
    }
  }

  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c(factor1, factor2))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- correction
  rownames(output$result) <- 1:dim(output$result)[1]
  class(output) <- "HRM"

  return (output)
}



#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor column name of the data frame X of within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between <- function(X, alpha, group , subgroup, factor, subject, data, testing = rep(1,7), formula, nonparametric, np.correction ){

  ranked <- NULL
  varQGlobal <- NULL
  correction <- NULL
  tmpQ1g <- NULL
  tmpQ2g <- NULL

  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  for(i in 1:7){
    if(testing[i]) {
      temp[[i]] <- hrm.2w.1f(X, alpha, group , subgroup, factor, subject, data, H = i, "", nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g )
    }
  }

  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group, subgroup), c(factor))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- correction
  rownames(output$result) <- 1:dim(output$result)[1]
  class(output) <- "HRM"

  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between.within <- function(X, alpha, group , subgroup, factor1, factor2, subject, data, testing = rep(1,15), formula, nonparametric, np.correction ){

  ranked <- NULL
  varQGlobal <- NULL
  correction <- NULL

  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  for(i in 1:15){
    if(testing[i]) {
      temp[[i]] <- hrm.2w.2f(X, alpha, group , subgroup, factor1, factor2, subject, data, H = i, "", nonparametric, ranked, varQGlobal, np.correction )
    }
  }
  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]], temp[[8]], temp[[9]], temp[[10]], temp[[11]], temp[[12]], temp[[13]], temp[[14]], temp[[15]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group, subgroup), c(factor1, factor2))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- correction
  rownames(output$result) <- 1:dim(output$result)[1]
  class(output) <- "HRM"
  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#'
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param factor3 column name of the data frame X of the third within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.3.between <- function(X, alpha, group , factor1, factor2, factor3, subject, data, testing = rep(1,15), formula, nonparametric, np.correction ){

  ranked <- NULL
  varQGlobal <- NULL
  correction <- NULL

  temp0 <- if(testing[1]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3, subject, data, "P", "J", "J", "J",  paste(as.character(group) ), nonparametric, ranked,varQGlobal, np.correction )}
  temp1 <- if(testing[2]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "P", "J", "J", paste(as.character(factor1)), nonparametric, ranked,varQGlobal, np.correction )}
  temp2 <- if(testing[3]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "J", "P", "J", paste(as.character(factor2)), nonparametric, ranked,varQGlobal, np.correction )}
  temp3 <- if(testing[4]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "J", "J", "P", paste(as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp4 <- if(testing[5]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "P", "J", "J", paste(as.character(group),":",as.character(factor1)), nonparametric, ranked,varQGlobal, np.correction )}
  temp5 <- if(testing[6]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "J", "P", "J", paste(as.character(group),":",as.character(factor2)), nonparametric, ranked,varQGlobal, np.correction )}
  temp6 <- if(testing[7]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "J", "J", "P", paste(as.character(group),":",as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp7 <- if(testing[8]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "P", "P", "J", paste(as.character(factor1),":",as.character(factor2)), nonparametric, ranked,varQGlobal, np.correction )}
  temp8 <- if(testing[9]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "P", "J", "P", paste(as.character(factor1),":",as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp9 <- if(testing[10]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "J", "P", "P",paste(as.character(factor2),":",as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp10 <- if(testing[11]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "P", "P", "J", paste(as.character(group),":",as.character(factor1), ":", as.character(factor2)), nonparametric, ranked,varQGlobal, np.correction )}
  temp11 <- if(testing[12]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "P", "J", "P", paste(as.character(group),":",as.character(factor1), ":", as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp12 <- if(testing[13]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "J", "P", "P", paste(as.character(group),":",as.character(factor2), ":", as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp13 <- if(testing[14]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "J", "P", "P", "P", paste(as.character(factor1),":",as.character(factor2), ":", as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}
  temp14 <- if(testing[15]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, "P", "P", "P", "P", paste(as.character(group),":",as.character(factor1), ":", as.character(factor2), ":", as.character(factor3)), nonparametric, ranked,varQGlobal, np.correction )}

  output <- list()
  output$result <- rbind(temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c(factor1, factor2, factor3))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  output$np.correction <- correction
  rownames(output$result) <- 1:dim(output$result)[1]
  class(output) <- "HRM"
  return (output)
}


#' Test for Multi-Factor High-Dimensional Repeated Measures
#'
#' @description Performing main and interaction effects of up to three whole- or subplot-factors. In total, a maximum of four factors can be used. There are two different S3 methods available. The first method requires a list of matrices in the wide table format. The second methodl requres a data.frame in the long table format.
#' @param data Either a data.frame (one observation per row) or a list with matrices (one subject per row) for all groups containing the data
#' @param formula A model formula object. The left hand side contains the response variable and the right hand side contains the whole- and subplot factors.
#' @param subject column name within the data frame X identifying the subjects
#' @param alpha alpha level used for calculating the critical value for the test
#' @return Returns an object from class HRM containing
#' @return \item{result}{A dataframe with the results from the hypotheses tests.}
#' @return \item{formula}{The formula object which was used.}
#' @return \item{alpha}{The type-I error rate which was used.}
#' @return \item{subject}{The column name identifying the subjects.}
#' @return \item{factors}{A list containing the whole- and subplot factors.}
#' @return \item{data}{The data.frame or list containing the data.}
#' @keywords internal
hrm_test_internal <- function(formula, data, alpha = 0.05,  subject, variable, nonparametric, np.correction){

  if(missing(data) || !is.data.frame(data)){
    stop("dataframe needed")
  }
  if(missing(subject) || !is.character(subject)){
    stop("subject column name not specified")
  }
  if(!is.double(alpha)){
    stop("alpha level needs to be a number between 0 and 1")
  }
  if(is.double(alpha)){
    if(alpha > 1 || alpha < 0){
      stop("alpha level needs to be a number between 0 and 1")
    }
  }

  # convert whole/subplot factor columns to type factor
  tryCatch({
    nfactors <- length(attr(terms.formula(formula), "variables"))
    for(i in 3:nfactors){
      cname <- as.character(attr(terms.formula(formula), "variables")[[i]])
      if( !is.factor( data[,cname] ) ) {
        data[,cname] <- as.factor(data[,cname])
      }
    }
  }, warning = function(w) "", error = function(e) { paste("One of the factor columns could not be converted to a factor variable." ) } )

  dat <- model.frame(formula, data)
  dat2 <- data.frame(dat,subj=data[,subject], variable_internal = data[, variable])

  m <- ncol(dat)

  if(!is.numeric(dat[,1])){
    stop("Response variable needs to be numeric!")
  }

  # find out, in which columns are the wholeplot or subplot factors
  s1<-subset(dat2, dat2$subj==dat2$subj[1])
  measurements <- dim(s1)[1]
  countSubplotFactor <- 1
  wholeplot<-rep(-1, m)
  subplot<-rep(-1, m)
  for(i in 2:m){
    if(!is.factor(dat2[,i])){
      stop(paste("The column ", colnames(dat2)[i], " is not a factor." ))
    }
    if(length(unique(s1[,i]))==nlevels(dat2[,i])){
      subplot[i]<-1
      countSubplotFactor <- countSubplotFactor*nlevels(s1[,i])
    }
    else{
      wholeplot[i]<-1
    }
  }
  wholeplot <- which(wholeplot==1)
  subplot <- which( subplot==1)

  if(is.null(variable)) {
    if(!(measurements == countSubplotFactor)){
      stop(paste("The number of repeated measurements per subject (", measurements, ") is uneqal to the number of levels of the subplot factors (", countSubplotFactor, ")."))
    }
  }
  if(!is.null(variable)) {
    p <- nlevels(dat2$variable_internal)
    if(!(measurements/p == countSubplotFactor)){
      stop(paste("The number of repeated measurements per subject (", measurements/p, ") is uneqal to the number of levels of the subplot factors (", countSubplotFactor, ")."))
    }
  }
  if(length(wholeplot)>2){
    stop("Too many factors are used! Only two wholelot-factors are supported.")
  }
  if(length(subplot)>6){
    stop("Too many factors are used! Only five subplotlot-factors are supported.")
  }
  if(length(wholeplot)>1 & length(subplot)==3){
    stop("Too many factors are used! Only one whole- and three subplot-factors are supported.")
  }
  if(length(subplot)<1 & length(wholeplot)>1){
    stop("The model needs at least one within-subject factor.")
  }
  if(length(subplot)>5 & length(wholeplot)<1){
    stop("The model supports up to five subplot-factor when using no wholeplot-factors.")
  }
  if((length(wholeplot) > 1  | length(subplot) > 1 ) & !is.null(variable)) {
    stop("The package currently supports Multivariate Repeated Measures only for models with 1 whole-plot and 1 sub-plot factor.")
  }

  # Case: no wholeplot, one subpot factor
  if(length(wholeplot) < 1 & length(subplot) == 1){
    factor1 <- colnames(dat2)[subplot[1]]
    x<-attributes(terms.formula(formula))$term.labels

    subplot<-colnames(dat2[,subplot])
    X<-data
    data <- colnames(dat)[1]
    if(is.null(variable)) {
      return(hrm.test.1.one(X, alpha , factor1, subject, data, formula, nonparametric, np.correction ))
    } else {
      return(hrm.mv.1w.1f(X, alpha, NULL , factor1, subject, data, variable, formula, nonparametric ))
    }
  }

  # Case: no wholeplot, two subpot factor
  if(length(wholeplot) < 1 & length(subplot) == 2){
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    x<-attributes(terms.formula(formula))$term.labels

    testing <- rep(0,3)
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 2 factors
      if(l == 2){testing[3]<-1}
      else{
        if(factor1 == x[i]){
          testing[1]<-1
        }
        else if(factor2 == x[i]){testing[2]<-1}
      }
    }

    subplot<-colnames(dat2[,subplot])
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.two(X, alpha , factor1, factor2, subject, data, formula, testing, nonparametric, np.correction ))

  }

  # Case: no wholeplot, three subpot factor
  if(length(wholeplot) < 1 & length(subplot) == 3){
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    factor3 <- colnames(dat2)[subplot[3]]
    x<-attributes(terms.formula(formula))$term.labels

    testing <- rep(0,2^3-1)
    testing[1:3] <- 1
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 2 factors
      if(l == 2){
        if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[4] <- 1
          } else {
            testing[5] <- 1
          }
        } else {
          testing[6] <- 1
        }
      }

      # interaction hypothesis of 3 factors
      if(l == 3){
        testing[7] <- 1
      }
    }

    subplot<-colnames(dat2[,subplot])
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.3.three(X, alpha , factor1, factor2, factor3, subject, data, formula, testing, nonparametric, np.correction ))

  }


  # Case: no wholeplot, four subpot factor
  if(length(wholeplot) < 1 & length(subplot) == 4){
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    factor3 <- colnames(dat2)[subplot[3]]
    factor4 <- colnames(dat2)[subplot[4]]
    x<-attributes(terms.formula(formula))$term.labels
    subplot<-colnames(dat2[,subplot])
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.4.four(X, alpha , factor1, factor2, factor3, factor4, subject, data, formula, testing = rep(1,2^4-1), nonparametric, np.correction ))

  }

  # Case: no wholeplot, five subpot factor
  if(length(wholeplot) < 1 & length(subplot) == 5){
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    factor3 <- colnames(dat2)[subplot[3]]
    factor4 <- colnames(dat2)[subplot[4]]
    factor5 <- colnames(dat2)[subplot[5]]
    x<-attributes(terms.formula(formula))$term.labels
    subplot<-colnames(dat2[,subplot])
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.5.five(X, alpha , factor1, factor2, factor3, factor4, factor5, subject, data, formula, testing = rep(1,2^5-1), nonparametric, np.correction ))

  }

  # Case: one wholeplot, no subpot factor
  if(length(wholeplot) == 1 & length(subplot) < 1){
    group <- colnames(dat2)[wholeplot[1]]
    x<-attributes(terms.formula(formula))$term.labels

    wholeplot<-colnames(dat2[,group])
    X<-data
    data <- colnames(dat)[1]
    if(is.null(variable)){
      return(hrm.test.1.none(X, alpha , group, subject, data, formula, nonparametric ))
    } else {
      return(hrm.mv.1w.1f(X, alpha, group , NULL, subject, data, variable, formula, nonparametric ))
    }
  }

  # Case: 1 whole and 1 subplot factor
  if(length(wholeplot)==1 & length(subplot)==1){
    group <- colnames(dat2)[wholeplot[1]]
    factor1 <- colnames(dat2)[subplot[1]]
    x<-attributes(terms.formula(formula))$term.labels

    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])

    testing <- rep(0,4)
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 2 factors
      if(l == 2){testing[4]<-1}

      else{
        if(group == x[i]){
          testing[1]<-1
          testing[2]<-1
        }
        else if(factor1 == x[i]){testing[3]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    if(is.null(variable)) {
      return(hrm.test.2.one(X, alpha, group , factor1, subject, data, testing, formula, nonparametric, np.correction ))
    } else {
      return(hrm.mv.1w.1f(X, alpha, group, factor1, subject, data, variable, formula, nonparametric))
    }
  }

  # Case: 2 wholeplot, 1 subplot factor
  if(length(wholeplot)==2 & length(subplot)==1){
    group <- colnames(dat2)[wholeplot[1]]
    subgroup <- colnames(dat2)[wholeplot[2]]
    factor1 <- colnames(dat2)[subplot[1]]
    x<-attributes(terms.formula(formula))$term.labels

    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])

    testing <- rep(0,7)
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 4 factors
      if(l == 3){testing[7]<-1}

      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            testing[4]<-1
          }
          if(grepl(factor1,x[i])){
            testing[5]<-1
          }
        }
        else if(grepl(subgroup,x[i])){
          if(grepl(factor1,x[i])){
            testing[6]<-1
          }
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]<-1}
        else if(subgroup == x[i]){testing[2]<-1}
        else if(factor1 == x[i]){testing[3]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.between(X, alpha, group , subgroup, factor1, subject, data, testing, formula, nonparametric, np.correction ))
  }

  # Case: 1 wholeplot, 2 subplot factors
  if(length(wholeplot)==1 & length(subplot)==2){
    group <- colnames(dat2)[wholeplot[1]]
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    x<-attributes(terms.formula(formula))$term.labels

    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])

    testing <- rep(0,7)
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 3 factors
      if(l == 3){testing[7]<-1}

      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            testing[4]<-1
          }
          if(grepl(factor2,x[i])){
            testing[5]<-1
          }
        }
        else if(grepl(factor1,x[i])){
          testing[6]<-1
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]<-1}
        else if(factor1 == x[i]){testing[2]<-1}
        else {testing[3]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.within(X, alpha, group , factor1, factor2, subject, data, testing, formula, nonparametric, np.correction ))
  }

  # Case: 1 wholeplot, 3 subplot factors
  if(length(wholeplot)==1 & length(subplot)==3){
    group <- colnames(dat2)[wholeplot[1]]
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    factor3 <- colnames(dat2)[subplot[3]]

    x<-attributes(terms.formula(formula))$term.labels

    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])

    testing <- rep(0,15)
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 3 factors
      if(l == 4){testing[15]<-1}

      # find out which interaction hypothesis of 3 factors is tested
      else if(l==3){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            if(grepl(factor2,x[i])){
              testing[11]<-1
            }
            else{
              testing[12]<-1
            }
          }
          if(grepl(factor2,x[i])){
            if(grepl(factor3,x[i])){
              testing[13]<-1
            }
          }
        }
        else if(grepl(factor1,x[i])){
          testing[14]<-1
        }
      }
      # l = 2
      else if (l==2){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            testing[5]<-1
          }
          if(grepl(factor2,x[i])){
            testing[6]<-1
          }
          if(grepl(factor3,x[i])){
            testing[7]<-1
          }
        }
        else if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[8]<-1
          }
          if(grepl(factor3,x[i])){
            testing[9]<-1
          }
        }
        else if(grepl(factor2,x[i])){
          if(grepl(factor3,x[i])){
            testing[10]<-1
          }
        }
      }
      else if(l==1){
        if(grepl(group, x[[i]])){
          testing[1]<-1
        }
        if(grepl(factor1, x[[i]])){
          testing[2]<-1
        }
        if(grepl(factor2, x[[i]])){
          testing[3]<-1
        }
        if(grepl(factor3, x[[i]])){
          testing[4]<-1
        }
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.3.between(X, alpha, group , factor1, factor2, factor3, subject, data, testing, formula, nonparametric, np.correction ))
  }



  if(length(wholeplot)==2 & length(subplot)==2){
    group <- colnames(dat2)[wholeplot[1]]
    subgroup <- colnames(dat2)[wholeplot[2]]
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    x<-attributes(terms.formula(formula))$term.labels

    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])

    testing <- rep(0,15)
    for(i in 1:length(x)){

      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])

      # interaction hypothesis of 4 factors
      if(l == 4){testing[15]<-1}

      # find out which interaction hypothesis of 3 factors is tested
      else if(l==3){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            if(grepl(factor1,x[i])){
              testing[11]<-1
            }
            else{
              testing[12]<-1
            }
          }
          else{
            testing[13]<-1
          }
        }
        else{
          testing[14]<-1
        }
      }

      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            testing[5]<-1
          }
          if(grepl(factor1,x[i])){
            testing[6]<-1
          }
          if(grepl(factor2,x[i])){
            testing[7]<-1
          }
        }
        else if(grepl(subgroup,x[i])){
          if(grepl(factor1,x[i])){
            testing[8]<-1
          }
          if(grepl(factor2,x[i])){
            testing[9]<-1
          }
        }
        else if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[10]<-1
          }
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]<-1}
        else if(subgroup == x[i]){testing[2]<-1}
        else if(factor1 == x[i]){testing[3]<-1}
        else {testing[4]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.between.within(X, alpha, group , subgroup, factor1, factor2, subject, data, testing, formula, nonparametric, np.correction ))
  }
}
