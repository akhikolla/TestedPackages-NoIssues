#' Create Historical Vertical Data Frames From Horizontal Data Frames
#'
#' \code{verticalize3()} returns a vertically formatted demographic data frame 
#' organized to create historical projection matrices, given a horizontally
#' formatted input data frame.
#'
#' @param data The horizontal data file.
#' @param noyears The number of years or observation periods in the dataset.
#' @param firstyear The first year or time of observation.
#' @param popidcol A variable name or column number corresponding to the identity
#' of the population for each individual.
#' @param patchidcol A variable name or column number corresponding to the identity
#' of the patch for each individual, if patches have been designated within
#' populations.
#' @param individcol A variable name or column number corresponding to the identity
#' of each individual.
#' @param blocksize The number of variables corresponding to each time in the
#' input dataset designated in \code{data}.
#' @param xcol A variable name or column number corresponding to the x coordinate
#' of each individual in Cartesian space.
#' @param ycol A variable name or column number corresponding to the y coordinate
#' of each individual in Cartesian space.
#' @param juvcol A variable name or column number that marks individuals in
#' immature stages within the dataset. The \code{verticalize3()} function assumes 
#' that immature individuals are identified in this variable marked with a number
#' equal to or greater than 1, and that mature individuals are marked as 0 or NA.
#' @param sizeacol A variable name or column number corresponding to the
#' size entry associated with the first year or observation time in the dataset.
#' @param sizebcol A second variable name or column number corresponding to the
#' size entry associated with the first year or observation time in the dataset.
#' @param sizeccol A third variable name or column number corresponding to the
#' size entry associated with the first year or observation time in the dataset.
#' @param repstracol A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, associated with the 
#' first year or observation period in the input dataset. This can be binomial or 
#' count data, and is used to in analysis of the probability of reproduction.
#' @param repstrbcol A second variable name or column number corresponding to
#' the production of reproductive structures, such as flowers, associated with
#' the first year or observation period in the input dataset. This can be binomial
#' or count data, and is used to in analysis of the probability of reproduction.
#' @param fecacol A variable name or column number denoting fecundity associated
#' with the first year or observation time in the input dataset. This may represent
#' egg counts, fruit counts, seed production, etc.
#' @param fecbcol A second variable name or column number denoting fecundity
#' associated with the first year or observation time in the input dataset. This 
#' may represent egg counts, fruit counts, seed production, etc.
#' @param indcovacol A variable name or column number corresponding to an
#' individual covariate to be used in analysis.
#' @param indcovbcol A variable name or column number corresponding to an
#' individual covariate to be used in analysis.
#' @param indcovccol A second variable name or column number corresponding to an
#' individual covariate to be used in analysis.
#' @param aliveacol A variable name or column number that provides information
#' on whether an individual is alive at a given time. If used, living status must
#' be designated as binomial (living = 1, dead = 0).
#' @param deadacol A variable name or column number that provides information on
#' whether an individual is alive at a given time. If used, dead status must be
#' designated as binomial (dead = 1, living = 0).
#' @param obsacol A variable name or column number providing information on whether
#' an individual is in an observable stage at a given time. If used, observation
#' status must be designated as binomial (observed = 1, not observed = 0).
#' @param nonobsacol A variable name or column number providing information on
#' whether an individual is in an unobservable stage at a given time. If used,
#' observation status must be designated as binomial (not observed = 1, 
#' observed = 0).
#' @param censorcol A variable name or column number corresponding to the first
#' entry of a censor variable, used to distinguish between entries to use and
#' entries not to use, or to designate entries with special issues that require
#' further attention. If used, this should be associated with the first year or
#' observation time, and all other years or times must also have censor columns.
#' @param repstrrel This is a scalar multiplier on variable \code{repstrbcol} to
#' make it equivalent to \code{repstracol}. This can be useful if two reproductive
#' status variables have related but unequal units, for example if \code{repstracol}
#' refers to one-flowered stems while \code{repstrbcol} refers to two-flowered
#' stems. Defaults to 1.
#' @param fecrel This is a scalar multiplier on variable \code{fecbcol} to make it
#' equivalent to \code{fecacol}. This can be useful if two fecundity variables have
#' related but unequal units. Defaults to 1.
#' @param stagecol Optional variable name or column number corresponding to life
#' history stage at a given time.
#' @param stageassign The stageframe object identifying the life history model
#' being operationalized. Note that if \code{stagecol} is provided, then this
#' stageframe is not used for stage designation.
#' @param stagesize A variable name or column number describing which size variable
#' to use in stage estimation. Defaults to NA, and can also take \code{sizea}, \code{sizeb},
#' \code{sizec}, or \code{sizeadded}, depending on which size variable is chosen.
#' @param censorkeep The value of the censor variable identifying data to be
#' included in analysis. Defaults to 0, but may take any value including NA.
#' @param censor A logical variable determining whether the output data should be 
#' censored using the variable defined in \code{censorcol}. Defaults to FALSE.
#' @param spacing The spacing at which density should be estimated, if density
#' estimation is desired and x and y coordinates are supplied. Given in the same
#' units as those used in the x and y coordinates given in \code{xcol} and \code{ycol}.
#' Defaults to NA.
#' @param NAas0 If TRUE, then all NA entries for size and fecundity variables will
#' be set to 0. This can help increase the sample size analyzed by \code{\link{modelsearch}()},
#' but should only be used when it is clear that this substitution is biologically
#' realistic. Defaults to FALSE.
#' @param NRasRep If TRUE, then will treat non-reproductive but mature individuals
#' as reproductive during stage assignment. This can be useful when a matrix is
#' desired without separation of reproductive and non-reproductive but mature
#' stages of the same size. Only used if \code{stageassign} is set to a stageframe.
#' Defaults to FALSE.
#' @param reduce A logical variable determining whether unused variables and some
#' invariant state variables should be removed from the output dataset. Defaults to
#' TRUE.
#' 
#' @return If all inputs are properly formatted, then this function will output a
#' historical vertical data frame (class \code{hfvdata}), meaning that the output
#' data frame will have three consecutive times of size and reproductive data per 
#' individual per row. This data frame is in standard format for all functions used
#' in \code{lefko3}, and so can be used without further modification.
#' 
#' @return Variables in this data frame include the following:
#' \item{rowid}{Unique identifier for the row of the data frame.}
#' \item{popid}{Unique identifier for the population, if given.}
#' \item{patchid}{Unique identifier for patch within population, if given.}
#' \item{individ}{Unique identifier for the individual.}
#' \item{year2}{Year or time at time \emph{t}.}
#' \item{firstseen}{Year or time of first observation.}
#' \item{lastseen}{Year or time of last observation.}
#' \item{obsage}{Observed age in time \emph{t}, assuming first observation
#' corresponds to age = 0.}
#' \item{obslifespan}{Observed lifespan, given as \code{lastseen - firstseen + 1}.}
#' \item{xpos1,xpos2,xpos3}{X position in Cartesian space in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively, if provided.}
#' \item{ypos1,ypos2,ypos3}{Y position in Cartesian space in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively, if provided.}
#' \item{sizea1,sizea2,sizea3}{Main size measurement in times \emph{t}-1, \emph{t}, and 
#' \emph{t}+1, respectively.}
#' \item{sizeb1,sizeb2,sizeb3}{Secondary size measurement in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively.}
#' \item{sizec1,sizec2,sizec3}{Tertiary measurement in times \emph{t}-1, \emph{t}, and
#' \emph{t}+1, respectively.}
#' \item{size1added,size2added,size3added}{Sum of primary, secondary, and tertiary
#' size measurements in timea \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstra1,repstra2,repstra3}{Main numbers of reproductive structures in
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstrb1,repstrb2,repstrb3}{Secondary numbers of reproductive structures
#' in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstr1added,repstr2added,repstr3added}{Sum of primary and secondary
#' reproductive structures in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{feca1,feca2,feca3}{Main numbers of offspring in times \emph{t}-1, \emph{t}, and 
#' \emph{t}+1, respectively.}
#' \item{fecb1,fecb2, fecb3}{Secondary numbers of offspring in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively.}
#' \item{fec1added,fec2added,fec3added}{Sum of primary and secondary fecundity in
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{censor1,censor2,censor3}{Censor state values in times \emph{t}-1, \emph{t}, and
#' \emph{t}+1, respectively.}
#' \item{juvgiven1,juvgiven2,juvgiven3}{Binomial variable indicating whether
#' individual is juvenile in times \emph{t}-1, \emph{t}, and \emph{t}+1. Only given if 
#' \code{juvcol} is provided.}
#' \item{obsstatus1,obsstatus2,obsstatus3}{Binomial observation state in times
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstatus1,repstatus2,repstatus3}{Binomial reproductive state in times 
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecstatus1,fecstatus2,fecstatus3}{Binomial offspring production state in
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{matstatus1,matstatus2,matstatus3}{Binomial maturity state in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{alive1,alive2,alive3}{Binomial state as alive in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively.}
#' \item{density}{Density of individuals per unit designated in \code{spacing}. Only
#' given if spacing is not NA.}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#'
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                          sizeacol = "Volume88", repstracol = "FCODE88",
#'                          fecacol = "Intactseed88", deadacol = "Dead1988",
#'                          nonobsacol = "Dormant1988", stageassign = lathframe,
#'                          stagesize = "sizea", censorcol = "Missing1988",
#'                          censorkeep = NA, censor = TRUE)
#' summary(lathvert)
#' 
#' @export
verticalize3 <- function(data, noyears, firstyear, popidcol = 0, patchidcol = 0, individcol= 0, 
                         blocksize, xcol = 0, ycol = 0, juvcol = 0, sizeacol, sizebcol = 0, 
                         sizeccol = 0, repstracol = 0, repstrbcol = 0, fecacol = 0, fecbcol = 0, 
                         indcovacol = 0, indcovbcol = 0, indcovccol = 0, aliveacol = 0, deadacol = 0, 
                         obsacol = 0, nonobsacol = 0, censorcol = 0, repstrrel = 1, fecrel = 1, 
                         stagecol = 0, stageassign = NA, stagesize = NA, censorkeep = 0, 
                         censor = FALSE, spacing = NA, NAas0 = FALSE, NRasRep = FALSE, reduce = TRUE) {
  
  stassign <- rowid <- alive2 <- indataset <- censor1 <- censor2 <- censor3 <- censbool <- NULL
  
  popid <- NA
  patchid <- NA
  individ <- NA
  
  #This first section tests the input for valid entries
  if (is.character(popidcol)) {
    if (is.element(popidcol, names(data))) {
      true.popidcol <- which(names(data) == popidcol)
      popidcol <- true.popidcol
    } else {stop("Please enter popidcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(patchidcol)) {
    if (is.element(patchidcol, names(data))) {
      true.patchidcol <- which(names(data) == patchidcol)
      patchidcol <- true.patchidcol
    } else {stop("Please enter patchidcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(individcol)) {
    if (is.element(individcol, names(data))) {
      true.individcol <- which(names(data) == individcol)
      individcol <- true.individcol
    } else {stop("Please enter individcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(xcol)) {
    if (is.element(xcol, names(data))) {
      true.xcol <- which(names(data) == xcol)
      xcol <- true.xcol
    } else {stop("Please enter xcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(ycol)) {
    if (is.element(ycol, names(data))) {
      true.ycol <- which(names(data) == ycol)
      ycol <- true.ycol
    } else {stop("Please enter ycol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(juvcol)) {
    if (is.element(juvcol, names(data))) {
      true.juvcol <- which(names(data) == juvcol)
      juvcol <- true.juvcol
    } else {stop("Please enter juvcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(sizeacol)) {
    if (is.element(sizeacol, names(data))) {
      true.sizeacol <- which(names(data) == sizeacol)
      sizeacol <- true.sizeacol
    } else {stop("Please enter sizeacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(sizebcol)) {
    if (is.element(sizebcol, names(data))) {
      true.sizebcol <- which(names(data) == sizebcol)
      sizebcol <- true.sizebcol
    } else {stop("Please enter sizebcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(sizeccol)) {
    if (is.element(sizeccol, names(data))) {
      true.sizeccol <- which(names(data) == sizeccol)
      sizeccol <- true.sizeccol
    } else {stop("Please enter sizeccol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(repstracol)) {
    if (is.element(repstracol, names(data))) {
      true.repstracol <- which(names(data) == repstracol)
      repstracol <- true.repstracol
    } else {stop("Please enter repstracol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(repstrbcol)) {
    if (is.element(repstrbcol, names(data))) {
      true.repstrbcol <- which(names(data) == repstrbcol)
      repstrbcol <- true.repstrbcol
    } else {stop("Please enter repstrbcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(fecacol)) {
    if (is.element(fecacol, names(data))) {
      true.fecacol <- which(names(data) == fecacol)
      fecacol <- true.fecacol
    } else {stop("Please enter fecacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(fecbcol)) {
    if (is.element(fecbcol, names(data))) {
      true.fecbcol <- which(names(data) == fecbcol)
      fecbcol <- true.fecbcol
    } else {stop("Please enter fecbcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(indcovacol)) {
    if (is.element(indcovacol, names(data))) {
      true.indcovacol <- which(names(data) == indcovacol)
      indcovacol <- true.indcovacol
    } else {stop("Please enter indcovacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(indcovbcol)) {
    if (is.element(indcovbcol, names(data))) {
      true.indcovbcol <- which(names(data) == indcovbcol)
      indcovbcol <- true.indcovbcol
    } else {stop("Please enter indcovbcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(indcovccol)) {
    if (is.element(indcovccol, names(data))) {
      true.indcovccol <- which(names(data) == indcovccol)
      indcovccol <- true.indcovccol
    } else {stop("Please enter indcovccol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(aliveacol)) {
    if (is.element(aliveacol, names(data))) {
      true.aliveacol <- which(names(data) == aliveacol)
      aliveacol <- true.aliveacol
    } else {stop("Please enter aliveacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(deadacol)) {
    if (is.element(deadacol, names(data))) {
      true.deadacol <- which(names(data) == deadacol)
      deadacol <- true.deadacol
    } else {stop("Please enter deadacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(obsacol)) {
    if (is.element(obsacol, names(data))) {
      true.obsacol <- which(names(data) == obsacol)
      obsacol <- true.obsacol
    } else {stop("Please enter obsacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(nonobsacol)) {
    if (is.element(nonobsacol, names(data))) {
      true.nonobsacol <- which(names(data) == nonobsacol)
      nonobsacol <- true.nonobsacol
    } else {stop("Please enter nonobsacol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(censorcol)) {
    if (is.element(censorcol, names(data))) {
      true.censorcol <- which(names(data) == censorcol)
      censorcol <- true.censorcol
    } else {stop("Please enter censorcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if(!all(is.na(stageassign))) {
    if(!is.element("stageframe", class(stageassign))) {
      stop("The stageassign option can only take NA or a stageframe object as input.", call. = FALSE)
    }
    if(length(intersect(stagesize, c("sizea", "sizeb", "sizec", "sizeadded"))) == 0) {
      stop("The stagesize option must equal NA, 'sizea', 'sizeb', 'sizec', or 'sizeadded'. No other values are permitted.", call. = FALSE)
    }
  }
  
  if (is.character(stagecol)) {
    if (is.element(stagecol, names(data))) {
      true.stagecol <- which(names(data) == stagecol)
      stagecol <- true.stagecol
    } else {stop("Please enter stagecol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  # Here we will modify our approach to verticalization based on the input stageframe
  if (!all(is.na(stageassign)) & stagecol == 0) {
    
    stassign <- TRUE 
    
    if (stagesize == "sizeadded") {
      stagesizecol <- 4
    } else if (stagesize == "sizec") {
      stagesizecol <-3
    } else if (stagesize == "sizeb") {
      stagesizecol <- 2
    } else {
      stagesizecol <- 1
    }
    
  } else {
    stassign <- FALSE
    
    stageassign <- as.data.frame(matrix(c(NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA), ncol = 16))
    names(stageassign) <- c("stagenames", "size", "repstatus", "obsstatus", "propstatus", "immstatus",
                            "matstatus", "indataset", "binhalfwidth_raw", "min_age", "max_age", "sizebin_min",
                            "sizebin_max", "sizebin_center", "sizebin_width", "comments")
    
    stagesizecol <- 0
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  if (censor) {
    if (!is.element(censorkeep, data[,censorcol])) {
      fullcenvec <- as.vector(apply(as.matrix(c(1:noyears)), 1, function(X) {c(data[,(censorcol + (X - 1) * blocksize)])}))
      
      if (!is.element(censorkeep, fullcenvec)) {
        stop("Please enter a valid value for censorkeep. This value should occur in the censor variable within the dataset.", 
             call. = FALSE)
      }
    }
  }
  
  if (is.na(censorkeep) & censor) { # This section checks to see if NA is the censor value to keep
    censbool <- TRUE
    censorkeep <- 0
  } else {
    censbool <- FALSE
  }
  
  popdatalist.new <- pfj(data, stageassign, noyears, firstyear, (popidcol - 1), (patchidcol - 1), (individcol - 1), 
                         blocksize, (xcol - 1), (ycol - 1), (juvcol - 1), (sizeacol - 1), (sizebcol - 1), 
                         (sizeccol - 1), (repstracol - 1), (repstrbcol - 1), (fecacol - 1), (fecbcol - 1),
                         (indcovacol - 1), (indcovbcol - 1), (indcovccol - 1), (aliveacol - 1), (deadacol - 1), 
                         (obsacol - 1), (nonobsacol - 1), (censorcol - 1), (stagecol - 1), repstrrel, 
                         fecrel, NAas0, NRasRep, stassign, stagesizecol, censbool)
  
  popdata <- do.call("cbind.data.frame", popdatalist.new)
  
  names(popdata) <- c("rowid", "popid", "patchid", "individ", "year2", "firstseen", "lastseen", "obsage", 
                      "obslifespan", "xpos1", "ypos1", "sizea1", "sizeb1", "sizec1", "size1added", "repstra1",
                      "repstrb1", "repstr1added", "feca1", "fecb1", "fec1added", "indcova1", "indcovb1", 
                      "indcovc1", "censor1", "juvgiven1", "obsstatus1", "repstatus1", "fecstatus1", 
                      "matstatus1", "alive1", "stage1", "stage1index", "xpos2", "ypos2", "sizea2", "sizeb2", 
                      "sizec2", "size2added", "repstra2", "repstrb2", "repstr2added", "feca2", "fecb2", 
                      "fec2added", "indcova2", "indcovb2", "indcovc2", "censor2", "juvgiven2", "obsstatus2", 
                      "repstatus2", "fecstatus2", "matstatus2", "alive2", "stage2", "stage2index", 
                      "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "size3added", "repstra3", "repstrb3", 
                      "repstr3added", "feca3", "fecb3", "fec3added", "indcova3", "indcovb3", "indcovc3", 
                      "censor3", "juvgiven3", "obsstatus3", "repstatus3", "fecstatus3", "matstatus3", 
                      "alive3", "stage3", "stage3index")
  
  rownames(popdata) <- c(1:dim(popdata)[1])
  
  popdatareal <- subset(popdata, subset = (alive2 == 1))
  
  if (any(!is.na(popdatareal$popid))) {popdatareal$popid <- as.factor(popdatareal$popid)}
  
  if (any(!is.na(popdatareal$patchid))) {popdatareal$patchid <- as.factor(popdatareal$patchid)}
  
  if (censor) {
    popdatareal <- subset(popdatareal, censor1 == censorkeep)
    popdatareal <- subset(popdatareal, censor2 == censorkeep)
    popdatareal <- subset(popdatareal, censor3 == censorkeep)
  }
  
  if (!is.na(spacing)) {
    popdatareal$density <- .density3(popdatareal, which(names(popdatareal) == "xpos2"),
                                     which(names(popdatareal) == "ypos2"),
                                     which(names(popdatareal) == "year2"), spacing)
  }
  
  if (reduce) {
    if (all(is.na(popdatareal$xpos1)) | length(unique(popdatareal$xpos1)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xpos1"))]}
    if (all(is.na(popdatareal$ypos1)) | length(unique(popdatareal$ypos1)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ypos1"))]}
    if (all(is.na(popdatareal$xpos2)) | length(unique(popdatareal$xpos2)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xpos2"))]}
    if (all(is.na(popdatareal$ypos2)) | length(unique(popdatareal$ypos2)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ypos2"))]}
    if (all(is.na(popdatareal$xpos3)) | length(unique(popdatareal$xpos3)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xpos3"))]}
    if (all(is.na(popdatareal$ypos3)) | length(unique(popdatareal$ypos3)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ypos3"))]}
    
    if (!is.na(censorkeep)) {
      if (censorcol > 0 & censor) {
        if (all(popdatareal$censor1 == popdatareal$censor1[1])) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor1"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(popdatareal$censor2 == popdatareal$censor2[1])) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor2"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(popdatareal$censor3 == popdatareal$censor3[1])) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor3"))]
        }
      }
    } else {
      if (censorcol > 0 & censor) {
        if (all(is.na(popdatareal$censor1))) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor1"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(is.na(popdatareal$censor2))) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor2"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(is.na(popdatareal$censor3))) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor3"))]
        }
      }
    }
    
    if (all(is.na(popdatareal$sizea1)) | length(unique(popdatareal$sizea1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizea1"))]
    }
    if (all(is.na(popdatareal$sizea2)) | length(unique(popdatareal$sizea2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizea2"))]
    }
    if (all(is.na(popdatareal$sizea3)) | length(unique(popdatareal$sizea3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizea3"))]
    }
    
    if (all(is.na(popdatareal$sizeb1)) | length(unique(popdatareal$sizeb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizeb1"))]
    }
    if (all(is.na(popdatareal$sizeb2)) | length(unique(popdatareal$sizeb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizeb2"))]
    }
    if (all(is.na(popdatareal$sizeb3)) | length(unique(popdatareal$sizeb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizeb3"))]
    }
    
    if (all(is.na(popdatareal$sizec1)) | length(unique(popdatareal$sizec1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizec1"))]
    }
    if (all(is.na(popdatareal$sizec2)) | length(unique(popdatareal$sizec2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizec2"))]
    }
    if (all(is.na(popdatareal$sizec3)) | length(unique(popdatareal$sizec3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizec3"))]
    }
    
    if (isTRUE(all.equal(popdatareal$size1added, popdatareal$sizea1))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size1added"))]
    } else if (all(is.na(popdatareal$size1added)) | length(unique(popdatareal$size1added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size1added"))]
    }
    if (isTRUE(all.equal(popdatareal$size2added, popdatareal$sizea2))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size2added"))]
    } else if (all(is.na(popdatareal$size2added)) | length(unique(popdatareal$size2added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size2added"))]
    }
    if (isTRUE(all.equal(popdatareal$size3added, popdatareal$sizea3))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size3added"))]
    } else if (all(is.na(popdatareal$size3added)) | length(unique(popdatareal$size3added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size3added"))]
    }
    
    if (all(is.na(popdatareal$repstra1)) | length(unique(popdatareal$repstra1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstra1"))]
    }
    if (all(is.na(popdatareal$repstra2)) | length(unique(popdatareal$repstra2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstra2"))]
    }
    if (all(is.na(popdatareal$repstra3)) | length(unique(popdatareal$repstra3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstra3"))]
    }
    
    if (all(is.na(popdatareal$repstrb1)) | length(unique(popdatareal$repstrb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstrb1"))]
    }
    if (all(is.na(popdatareal$repstrb2)) | length(unique(popdatareal$repstrb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstrb2"))]
    }
    if (all(is.na(popdatareal$repstrb3)) | length(unique(popdatareal$repstrb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstrb3"))]
    }
    
    if (isTRUE(all.equal(popdatareal$repstr1added, popdatareal$repstr1a))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr1added"))]
    } else if (all(is.na(popdatareal$repstr1added)) | length(unique(popdatareal$repstr1added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr1added"))]
    }
    if (isTRUE(all.equal(popdatareal$repstr2added, popdatareal$repstr2a))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr2added"))]
    } else if (all(is.na(popdatareal$repstr2added)) | length(unique(popdatareal$repstr2added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr2added"))]
    }
    if (isTRUE(all.equal(popdatareal$repstr3added, popdatareal$repstr3a))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr3added"))]
    } else if (all(is.na(popdatareal$repstr3added)) | length(unique(popdatareal$repstr3added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr3added"))]
    }
    
    if (all(is.na(popdatareal$feca1)) | length(unique(popdatareal$feca1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="feca1"))]
    }
    if (all(is.na(popdatareal$feca2)) | length(unique(popdatareal$feca2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="feca2"))]
    }
    if (all(is.na(popdatareal$feca3)) | length(unique(popdatareal$feca3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="feca3"))]
    }
    
    if (all(is.na(popdatareal$fecb1)) | length(unique(popdatareal$fecb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecb1"))]
    }
    if (all(is.na(popdatareal$fecb2)) | length(unique(popdatareal$fecb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecb2"))]
    }
    if (all(is.na(popdatareal$fecb3)) | length(unique(popdatareal$fecb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecb3"))]
    }
    
    if (isTRUE(all.equal(popdatareal$fec1added, popdatareal$feca1))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec1added"))]
    } else if (all(is.na(popdatareal$fec1added)) | length(unique(popdatareal$fec1added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec1added"))]
    }
    if (isTRUE(all.equal(popdatareal$fec2added, popdatareal$feca2))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec2added"))]
    } else if (all(is.na(popdatareal$fec2added)) | length(unique(popdatareal$fec2added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec2added"))]
    }
    if (isTRUE(all.equal(popdatareal$fec3added, popdatareal$feca3))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec3added"))]
    } else if (all(is.na(popdatareal$fec3added)) | length(unique(popdatareal$fec3added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec3added"))]
    }
    
    if (all(is.na(popdatareal$indcova1)) | length(unique(popdatareal$indcova1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcova1"))]
    }
    if (all(is.na(popdatareal$indcova2)) | length(unique(popdatareal$indcova2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcova2"))]
    }
    if (all(is.na(popdatareal$indcova3)) | length(unique(popdatareal$indcova3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcova3"))]
    }
    
    if (all(is.na(popdatareal$indcovb1)) | length(unique(popdatareal$indcovb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcovb1"))]
    }
    if (all(is.na(popdatareal$indcovb2)) | length(unique(popdatareal$indcovb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcovb2"))]
    }
    if (all(is.na(popdatareal$indcovb3)) | length(unique(popdatareal$indcovb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcovb3"))]
    }
    
    if (all(is.na(popdatareal$indcovc1)) | length(unique(popdatareal$indcovc1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcovc1"))]
    }
    if (all(is.na(popdatareal$indcovc2)) | length(unique(popdatareal$indcovc2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcovc2"))]
    }
    if (all(is.na(popdatareal$indcovc3)) | length(unique(popdatareal$indcovc3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="indcovc3"))]
    }
    
    if (all(popdatareal$obsstatus1 == popdatareal$obsstatus1[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="obsstatus1"))]
    }
    if (all(popdatareal$obsstatus2 == popdatareal$obsstatus2[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="obsstatus2"))]
    }
    if (all(popdatareal$obsstatus3 == popdatareal$obsstatus3[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="obsstatus3"))]
    }
    
    if (all(popdatareal$repstatus1 == popdatareal$repstatus1[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstatus1"))]
    }
    if (all(popdatareal$repstatus2 == popdatareal$repstatus2[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstatus2"))]
    }
    if (all(popdatareal$repstatus3 == popdatareal$repstatus3[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstatus3"))]
    }
    
    if (all(popdatareal$fecstatus1 == popdatareal$fecstatus1[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecstatus1"))]
    }
    if (all(popdatareal$fecstatus2 == popdatareal$fecstatus2[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecstatus2"))]
    }
    if (all(popdatareal$fecstatus3 == popdatareal$fecstatus3[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecstatus3"))]
    }
    
    if (!censor) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor1"))]
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor2"))]
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor3"))]
    }
  }
  
  class(popdatareal) <- append(class(popdatareal), "hfvdata")
  
  return(popdatareal)
}

#' Create Historical Vertical Data Frame From Ahistorical Vertical Data Frame
#' 
#' \code{historicalize3()} returns a vertically formatted demographic data frame
#' organized to create historical projection matrices, given a vertically but
#' ahistorically formatted data frame. This data frame is in standard \code{lefko3}
#' format and can be used in all functions in the package.
#'
#' @param data The horizontal data file.
#' @param popidcol A variable name or column number corresponding to the identity
#' of the population for each individual.
#' @param patchidcol A variable name or column number corresponding to the identity
#' of the patch for each individual, if patches have been designated within
#' populations.
#' @param individcol A variable name or column number corresponding to the unique
#' identity of each individual.
#' @param year2col A variable name or column number corresponding to the year or
#' time in time \emph{t}.
#' @param year3col A variable name or column number corresponding to the year or
#' time in time \emph{t}+1.
#' @param xcol A variable name or column number corresponding to the x coordinate
#' of each individual in Cartesian space.
#' @param ycol A variable name or column number corresponding to the y coordinate
#' of each individual in Cartesian space.
#' @param sizea2col A variable name or column number corresponding to the primary
#' size entry in time \emph{t}.
#' @param sizea3col A variable name or column number corresponding to the primary 
#' size entry in time \emph{t}+1.
#' @param sizeb2col A variable name or column number corresponding to the
#' secondary size entry in time \emph{t}.
#' @param sizeb3col A variable name or column number corresponding to the
#' secondary size entry in time \emph{t}+1.
#' @param sizec2col A variable name or column number corresponding to the tertiary
#' size entry in time \emph{t}.
#' @param sizec3col A variable name or column number corresponding to the tertiary
#' size entry in time \emph{t}+1.
#' @param repstra2col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in time \emph{t}. This can
#' be binomial or count data, and is used to in analysis of the probability of 
#' reproduction.
#' @param repstra3col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in time \emph{t}+1. This
#' can be binomial or count data, and is used to in analysis of the probability of
#' reproduction.
#' @param repstrb2col A second variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in time \emph{t}. This can
#' be binomial or count data.
#' @param repstrb3col A second variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in time \emph{t}+1. This can
#' be binomial or count data.
#' @param feca2col A variable name or column number corresponding to fecundity in
#' time \emph{t}. This may represent egg counts, fruit counts, seed production,
#' etc.
#' @param feca3col A variable name or column number corresponding to fecundity in
#' time \emph{t}+1. This may represent egg counts, fruit counts, seed production,
#' etc.
#' @param fecb2col A second variable name or column number corresponding to 
#' fecundity in time \emph{t}. This may represent egg counts, fruit counts, seed
#' production, etc.
#' @param fecb3col A second variable name or column number corresponding to 
#' fecundity in time \emph{t}+1. This may represent egg counts, fruit counts, seed
#' production, etc.
#' @param indcova2col A variable name or column number corresponding to an
#' individual covariate to be used in analysis, in time \emph{t}.
#' @param indcova3col A variable name or column number corresponding to an
#' individual covariate to be used in analysis, in time \emph{t}+1.
#' @param indcovb2col A second variable name or column number corresponding to an
#' individual covariate to be used in analysis, in time \emph{t}.
#' @param indcovb3col A second variable name or column number corresponding to an
#' individual covariate to be used in analysis, in time \emph{t}+1.
#' @param indcovc2col A third variable name or column number corresponding to an
#' individual covariate to be used in analysis, in time \emph{t}.
#' @param indcovc3col A third variable name or column number corresponding to an
#' individual covariate to be used in analysis, in time \emph{t}+1.
#' @param alive2col A variable name or column number that provides information on
#' whether an individual is alive in time \emph{t}. If used, living status must be
#' designated as binomial (living = 1, dead = 0).
#' @param alive3col A variable name or column number that provides information on
#' whether an individual is alive in time \emph{t}+1. If used, living status must
#' be designated as binomial (living = 1, dead = 0).
#' @param dead2col A variable name or column number that provides information on
#' whether an individual is dead in time \emph{t}. If used, dead status
#' must be designated as binomial (dead = 1, living = 0).
#' @param dead3col A variable name or column number that provides information on
#' whether an individual is dead in time \emph{t}+1. If used, dead status must be
#' designated as binomial (dead = 1, living = 0).
#' @param obs2col A variable name or column number providing information on whether
#' an individual is in an observable stage in time \emph{t}. If used, observation
#' status must be designated as binomial (observed = 1, not observed = 0).
#' @param obs3col A variable name or column number providing information on whether
#' an individual is in an observable stage in time \emph{t}+1. If used, observation
#' status must be designated as binomial (observed = 1, not observed = 0).
#' @param nonobs2col A variable name or column number providing information on
#' whether an individual is in an unobservable stage in time \emph{t}. If used,
#' observation status must be designated as binomial (not observed = 1,
#' observed = 0).
#' @param nonobs3col A variable name or column number providing information on
#' whether an individual is in an unobservable stage in time \emph{t}+1. If used,
#' observation status must be designated as binomial (not observed = 1,
#' observed = 0).
#' @param repstrrel This is a scalar multiplier to make the variable represented by
#' \code{repstrb2col} equivalent to the variable represented by \code{repstra2col}.
#' This can be useful if two reproductive status variables have related but unequal
#' units, for example if \code{repstrb2col} refers to one-flowered stems while 
#' \code{repstra2col} refers to two-flowered stems.
#' @param fecrel This is a scalar multiplier that makes the variable represented by
#' \code{fecb2col} equivalent to the variable represented by \code{feca2col}. This can
#' be useful if two fecundity variables have related but unequal units.
#' @param stage2col Optional variable name or column number corresponding to life
#' history stage in time \emph{t}.
#' @param stage3col Optional variable name or column number corresponding to life
#' history stage in time \emph{t}+1.
#' @param juv2col A variable name or column number that marks individuals in
#' immature stages in time \emph{t}. The \code{historicalize3()} function assumes 
#' that immature individuals are identified in this variable marked with a number
#' equal to or greater than 1, and that mature individuals are marked as 0 or NA.
#' @param juv3col A variable name or column number that marks individuals in
#' immature stages in time \emph{t}+1. The \code{historicalize3()} function assumes 
#' that immature individuals are identified in this variable marked with a number
#' equal to or greater than 1, and that mature individuals are marked as 0 or NA.
#' @param stageassign The stageframe object identifying the life history model
#' being operationalized. Note that if \code{stage2col} is provided, then this
#' stageframe is not utilized in stage designation.
#' @param stagesize A variable name or column number describing which size variable
#' to use in stage estimation. Defaults to NA, and can also take \code{sizea}, \code{sizeb},
#' \code{sizec}, or \code{sizeadded}, depending on which size variable is chosen.
#' @param censor A logical variable determining whether the output data should be
#' censored using the variable defined in \code{censorcol}. Defaults to FALSE.
#' @param censorcol A variable name or column number corresponding to a censor
#' variable within the dataset, used to distinguish between entries to use and
#' those to discard from analysis, or to designate entries with special issues that
#' require further attention.
#' @param censorkeep The value of the censoring variable identifying data that
#' should be included in analysis. Defaults to 0, but may take any value including
#' NA.
#' @param spacing The spacing at which density should be estimated, if density
#' estimation is desired and x and y coordinates are supplied. Given in the same
#' units as those used in the x and y coordinates given in \code{xcol} and \code{ycol}.
#' Defaults to NA.
#' @param NAas0 If TRUE, then all NA entries for size and fecundity variables will
#' be set to 0. This can help increase the sample size analyzed by \code{\link{modelsearch}()},
#' but should only be used when it is clear that this substitution is biologically
#' realistic. Defaults to FALSE.
#' @param NRasRep If set to TRUE, then this function will treat non-reproductive
#' but mature individuals as reproductive during stage zssignment. This can be 
#' useful when a matrix is desired without separation of reproductive and
#' non-reproductive but mature stages of the same size. Only used if \code{stageassign}
#' is set to a stageframe. Defaults to FALSE.
#' @param reduce A logical variable determining whether unused variables and some
#' invariant state variables should be removed from the output dataset. Defaults
#' to TRUE.
#'
#' @return If all inputs are properly formatted, then this function will output a
#' historical vertical data frame (class \code{hfvdata}), meaning that the output
#' data frame will have three consecutive years of size and reproductive data per
#' individual per row. This data frame is in standard format for all functions used
#' in \code{lefko3}, and so can be used without further modification. Note that
#' determination of state in times *t*-1 and *t*+1 gives preference to condition in
#' time *t* within the input dataset. Conflicts in condition in input datasets that
#' have both times *t* and *t*+1 listed per row are resolved by using condition in
#' time *t*.
#' \item{rowid}{Unique identifier for the row of the data frame.}
#' \item{popid}{Unique identifier for the population, if given.}
#' \item{patchid}{Unique identifier for patch within population, if given.}
#' \item{individ}{Unique identifier for the individual.}
#' \item{year2}{Year or time at time \emph{t}.}
#' \item{firstseen}{Year or time of first observation.}
#' \item{lastseen}{Year or time of last observation.}
#' \item{obsage}{Observed age in time \emph{t}, assuming first observation
#' corresponds to age = 0.}
#' \item{obslifespan}{Observed lifespan, given as \code{lastseen - firstseen + 1}.}
#' \item{xpos1,xpos2,xpos3}{X position in Cartesian space in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively, if provided.}
#' \item{ypos1,ypos2,ypos3}{Y position in Cartesian space in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively, if provided.}
#' \item{sizea1,sizea2,sizea3}{Main size measurement in times \emph{t}-1, \emph{t}, and
#' \emph{t}+1, respectively.}
#' \item{sizeb1,sizeb2,sizeb3}{Secondary size measurement in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively.}
#' \item{sizec1,sizec2,sizec3}{Tertiary measurement in times \emph{t}-1, \emph{t}, and 
#' \emph{t}+1, respectively.}
#' \item{size1added,size2added,size3added}{Sum of primary, secondary, and tertiary
#' size measurements in timea \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstra1,repstra2,repstra3}{Main numbers of reproductive structures in
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstrb1,repstrb2,repstrb3}{Secondary numbers of reproductive structures
#' in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstr1added,repstr2added,repstr3added}{Sum of primary and secondary
#' reproductive structures in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{feca1,feca2,feca3}{Main numbers of offspring in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively.}
#' \item{fecb1,fecb2, fecb3}{Secondary numbers of offspring in times \emph{t}-1,
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{fec1added,fec2added,fec3added}{Sum of primary and secondary fecundity in
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{censor1,censor2,censor3}{Censor state values in times \emph{t}-1, \emph{t}, and
#' \emph{t}+1, respectively.}
#' \item{juvgiven1,juvgiven2,juvgiven3}{Binomial variable indicating whether
#' individual is juvenile in times \emph{t}-1, \emph{t}, and \emph{t}+1. Only given
#' if \code{juvcol} is provided.}
#' \item{obsstatus1,obsstatus2,obsstatus3}{Binomial observation state in times \emph{t}-1,
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstatus1,repstatus2,repstatus3}{Binomial reproductive state in times
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecstatus1,fecstatus2,fecstatus3}{Binomial offspring production state in
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{matstatus1,matstatus2,matstatus3}{Binomial maturity state in times \emph{t}-1,
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{alive1,alive2,alive3}{Binomial state as alive in times \emph{t}-1, \emph{t},
#' and \emph{t}+1, respectively.}
#' \item{density}{Density of individuals per unit designated in \code{spacing}. Only
#' given if spacing is not NA.}
#' 
#' @examples
#' data(cypvert)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg")
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
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", individcol = "plantid",
#'                             year2col = "year2", sizea2col = "Inf2.2", sizea3col = "Inf2.3",
#'                             sizeb2col = "Inf.2", sizeb3col = "Inf.3", sizec2col = "Veg.2",
#'                             sizec3col = "Veg.3", repstra2col = "Inf2.2", repstra3col = "Inf2.3",
#'                             repstrb2col = "Inf.2", repstrb3col = "Inf.3", feca2col = "Pod.2",
#'                             feca3col = "Pod.3", repstrrel = 2, stageassign = cypframe_raw,
#'                             stagesize = "sizeadded", censorcol = "censor", censor = FALSE,
#'                             NAas0 = TRUE, NRasRep = TRUE)
#' summary(cypraw_v2)
#'
#' @export
historicalize3 <- function(data, popidcol = 0, patchidcol = 0, individcol, year2col = 0,
                           year3col = 0, xcol = 0, ycol = 0, sizea2col = 0, sizea3col = 0,
                           sizeb2col = 0, sizeb3col = 0, sizec2col = 0, sizec3col = 0,
                           repstra2col = 0, repstra3col = 0, repstrb2col = 0,
                           repstrb3col = 0, feca2col = 0, feca3col = 0,  fecb2col = 0,
                           fecb3col = 0, indcova2col = 0, indcova3col = 0, indcovb2col = 0,
                           indcovb3col = 0, indcovc2col = 0, indcovc3col = 0,
                           alive2col = 0, alive3col = 0, dead2col = 0,
                           dead3col = 0, obs2col = 0, obs3col = 0, nonobs2col = 0,
                           nonobs3col = 0, repstrrel = 1, fecrel = 1, stage2col = 0,
                           stage3col = 0, juv2col = 0, juv3col = 0, stageassign = NA,
                           stagesize = NA, censor = FALSE, censorcol = 0, censorkeep = 0,
                           spacing = NA, NAas0 = FALSE, NRasRep = FALSE, reduce = TRUE) {
  
  alive2 <- indataset <- censor1 <- censor2 <- censor3 <- censbool <- NULL
  
  if (is.na(individcol)) {
    stop("Individual ID variable is required.", .call = FALSE)
  }
  
  if (is.na(year2col) & is.na(year3col)) {
    stop("Variable identifying either year2 (time t) or year3 (time t+1) is required.", .call = FALSE)
  }
  
  if (is.character(popidcol)) {
    if (is.element(popidcol, names(data))) {
      true.popidcol <- which(names(data) == popidcol)
      popidcol <- true.popidcol
    } else {
      stop("Please enter popidcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(patchidcol)) {
    if (is.element(patchidcol, names(data))) {
      true.patchidcol <- which(names(data) == patchidcol)
      patchidcol <- true.patchidcol
    } else {
      stop("Please enter patchidcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(individcol)) {
    if (is.element(individcol, names(data))) {
      true.individcol <- which(names(data) == individcol)
      individcol <- true.individcol
    } else {
      stop("Please enter individcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (year2col != 0) {
    if (is.character(year2col)) {
      if (is.element(year2col, names(data))) {
        true.year2col <- which(names(data) == year2col)
        year2col <- true.year2col
      } else {
        stop("Please enter year2col exactly as it appears in the dataset.", call. = FALSE)
      }
    }
  } 
  
  if (year3col != 0) {
    if (is.character(year3col)) {
      if (is.element(year3col, names(data))) {
        true.year3col <- which(names(data) == year3col)
        year3col <- true.year3col
      } else {
        stop("Please enter year3col exactly as it appears in the dataset.", call. = FALSE)
      }
    }
  }
  
  if (year2col != 0 & year3col == 0) {
    data$year3 <- data[,year2col] + 1
    year3col <- which(names(data) == "year3")
  } else if (year2col == 0 & year3col != 0) {
    data$year2 <- data[,year3col] - 1
    year2col <- which(names(data) == "year2")
  }
  
  if (is.character(xcol)) {
    if (is.element(xcol, names(data))) {
      true.xcol <- which(names(data) == xcol)
      xcol <- true.xcol
    } else {
      stop("Please enter xcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(ycol)) {
    if (is.element(ycol, names(data))) {
      true.ycol <- which(names(data) == ycol)
      ycol <- true.ycol
    } else {
      stop("Please enter ycol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizea2col)) {
    if (is.element(sizea2col, names(data))) {
      true.sizea2col <- which(names(data) == sizea2col)
      sizea2col <- true.sizea2col
    } else {
      stop("Please enter sizea2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizea3col)) {
    if (is.element(sizea3col, names(data))) {
      true.sizea3col <- which(names(data) == sizea3col)
      sizea3col <- true.sizea3col
    } else {
      stop("Please enter sizea3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizeb2col)) {
    if (is.element(sizeb2col, names(data))) {
      true.sizeb2col <- which(names(data) == sizeb2col)
      sizeb2col <- true.sizeb2col
    } else {
      stop("Please enter sizeb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizeb3col)) {
    if (is.element(sizeb3col, names(data))) {
      true.sizeb3col <- which(names(data) == sizeb3col)
      sizeb3col <- true.sizeb3col
    } else {
      stop("Please enter sizeb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizec2col)) {
    if (is.element(sizec2col, names(data))) {
      true.sizec2col <- which(names(data) == sizec2col)
      sizec2col <- true.sizec2col
    } else {
      stop("Please enter sizec2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizec3col)) {
    if (is.element(sizec3col, names(data))) {
      true.sizec3col <- which(names(data) == sizec3col)
      sizec3col <- true.sizec3col
    } else {
      stop("Please enter sizec3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstra2col)) {
    if (is.element(repstra2col, names(data))) {
      true.repstra2col <- which(names(data) == repstra2col)
      repstra2col <- true.repstra2col
    } else {
      stop("Please enter repstra2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstra3col)) {
    if (is.element(repstra3col, names(data))) {
      true.repstra3col <- which(names(data) == repstra3col)
      repstra3col <- true.repstra3col
    } else {
      stop("Please enter repstra3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstrb2col)) {
    if (is.element(repstrb2col, names(data))) {
      true.repstrb2col <- which(names(data) == repstrb2col)
      repstrb2col <- true.repstrb2col
    } else {
      stop("Please enter repstrb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstrb3col)) {
    if (is.element(repstrb3col, names(data))) {
      true.repstrb3col <- which(names(data) == repstrb3col)
      repstrb3col <- true.repstrb3col
    } else {
      stop("Please enter repstrb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(feca2col)) {
    if (is.element(feca2col, names(data))) {
      true.feca2col <- which(names(data) == feca2col)
      feca2col <- true.feca2col
    } else {
      stop("Please enter feca2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(feca3col)) {
    if (is.element(feca3col, names(data))) {
      true.feca3col <- which(names(data) == feca3col)
      feca3col <- true.feca3col
    } else {
      stop("Please enter feca3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(fecb2col)) {
    if (is.element(fecb2col, names(data))) {
      true.fecb2col <- which(names(data) == fecb2col)
      fecb2col <- true.fecb2col
    } else {
      stop("Please enter fecb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(fecb3col)) {
    if (is.element(fecb3col, names(data))) {
      true.fecb3col <- which(names(data) == fecb3col)
      fecb3col <- true.fecb3col
    } else {
      stop("Please enter fecb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(indcova2col)) {
    if (is.element(indcova2col, names(data))) {
      true.indcova2col <- which(names(data) == indcova2col)
      indcova2col <- true.indcova2col
    } else {
      stop("Please enter indcova2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(indcova3col)) {
    if (is.element(indcova3col, names(data))) {
      true.indcova3col <- which(names(data) == indcova3col)
      indcova3col <- true.indcova3col
    } else {
      stop("Please enter indcova3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(indcovb2col)) {
    if (is.element(indcovb2col, names(data))) {
      true.indcovb2col <- which(names(data) == indcovb2col)
      indcovb2col <- true.indcovb2col
    } else {
      stop("Please enter indcovb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(indcovb3col)) {
    if (is.element(indcovb3col, names(data))) {
      true.indcovb3col <- which(names(data) == indcovb3col)
      indcovb3col <- true.indcovb3col
    } else {
      stop("Please enter indcovb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(indcovc2col)) {
    if (is.element(indcovc2col, names(data))) {
      true.indcovc2col <- which(names(data) == indcovc2col)
      indcovc2col <- true.indcovc2col
    } else {
      stop("Please enter indcovc2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(indcovc3col)) {
    if (is.element(indcovc3col, names(data))) {
      true.indcovc3col <- which(names(data) == indcovc3col)
      indcovc3col <- true.indcovc3col
    } else {
      stop("Please enter indcovc3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(alive2col)) {
    if (is.element(alive2col, names(data))) {
      true.alive2col <- which(names(data) == alive2col)
      alive2col <- true.alive2col
    } else {
      stop("Please enter alive2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(alive3col)) {
    if (is.element(alive3col, names(data))) {
      true.alive3col <- which(names(data) == alive3col)
      alive3col <- true.alive3col
    } else {
      stop("Please enter alive3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(dead2col)) {
    if (is.element(dead2col, names(data))) {
      true.dead2col <- which(names(data) == dead2col)
      dead2col <- true.dead2col
    } else {
      stop("Please enter dead2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(dead3col)) {
    if (is.element(dead3col, names(data))) {
      true.dead3col <- which(names(data) == dead3col)
      dead3col <- true.dead3col
    } else {
      stop("Please enter dead3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(obs2col)) {
    if (is.element(obs2col, names(data))) {
      true.obs2col <- which(names(data) == obs2col)
      obs2col <- true.obs2col
    } else {
      stop("Please enter obs2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(obs3col)) {
    if (is.element(obs3col, names(data))) {
      true.obs3col <- which(names(data) == obs3col)
      obs3col <- true.obs3col
    } else {
      stop("Please enter obs3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(nonobs2col)) {
    if (is.element(nonobs2col, names(data))) {
      true.nonobs2col <- which(names(data) == nonobs2col)
      nonobs2col <- true.nonobs2col
    } else {
      stop("Please enter nonobs2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(nonobs3col)) {
    if (is.element(nonobs3col, names(data))) {
      true.nonobs3col <- which(names(data) == nonobs3col)
      nonobs3col <- true.nonobs3col
    } else {
      stop("Please enter nonobs3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(stage2col)) {
    if (is.element(stage2col, names(data))) {
      true.stage2col <- which(names(data) == stage2col)
      stage2col <- true.stage2col
    } else {
      stop("Please enter stage2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(stage3col)) {
    if (is.element(stage3col, names(data))) {
      true.stage3col <- which(names(data) == stage3col)
      stage3col <- true.stage3col
    } else {
      stop("Please enter stage3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(juv2col)) {
    if (is.element(juv2col, names(data))) {
      true.juv2col <- which(names(data) == juv2col)
      juv2col <- true.juv2col
    } else {
      stop("Please enter juv2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(juv3col)) {
    if (is.element(juv3col, names(data))) {
      true.juv3col <- which(names(data) == juv3col)
      juv3col <- true.juv3col
    } else {
      stop("Please enter juv3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(censorcol)) {
    if (is.element(censorcol, names(data))) {
      true.censorcol <- which(names(data) == censorcol)
      censorcol <- true.censorcol
    } else {
      stop("Please enter censorcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  if (!all(is.na(stageassign))) {
    
    stassign <- TRUE 
    
    if (stagesize == "sizeadded") {
      stagesizecol <- 4
    } else if (stagesize == "sizec") {
      stagesizecol <-3
    } else if (stagesize == "sizeb") {
      stagesizecol <- 2
    } else {
      stagesizecol <- 1
    }
    
  } else {
    stassign <- FALSE
    
    stageassign <- as.data.frame(matrix(c(NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA), ncol = 16))
    names(stageassign) <- c("stagenames", "size", "repstatus", "obsstatus", "propstatus", "immstatus",
                            "matstatus", "indataset", "binhalfwidth_raw", "min_age", "max_age", "sizebin_min",
                            "sizebin_max", "sizebin_center", "sizebin_width", "comments")
    
    stagesizecol <- 0
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  if (censor) {
    if (!is.element(censorkeep, data[,censorcol])) {
      stop("Please enter a valid value for censorkeep. This value should occur in the censor variable within the dataset.", 
           call. = FALSE)
    }
  }
  
  if (is.na(censorkeep) & censor) {
    censbool <- TRUE
    censorkeep <- 0
  } else {
    censbool <- FALSE
  }
  
  
  
  
  
  
  
  
  
  
  
  
  popdatalist.new <- jpf(data, stageassign, (popidcol - 1), (patchidcol - 1), (individcol - 1),
                         (year2col - 1), (year3col - 1), (xcol - 1), (ycol - 1), (juv2col - 1), (juv3col - 1),
                         (sizea2col - 1), (sizea3col - 1), (sizeb2col - 1), (sizeb3col - 1), (sizec2col - 1),
                         (sizec3col - 1), (repstra2col - 1), (repstra3col - 1), (repstrb2col - 1), 
                         (repstrb3col - 1), (feca2col - 1), (feca3col - 1), (fecb2col - 1), (fecb3col - 1), 
                         (indcova2col - 1), (indcova3col - 1), (indcovb2col - 1), (indcovb3col - 1), 
                         (indcovc2col - 1), (indcovc3col - 1), (alive2col - 1), (alive3col - 1), 
                         (dead2col - 1), (dead3col - 1), (obs2col - 1), (obs3col - 1), (nonobs2col - 1), 
                         (nonobs3col - 1), repstrrel, fecrel, (stage2col - 1), (stage3col - 1), 
                         (censorcol - 1), NAas0, NRasRep, stassign, stagesizecol, censbool)
  
  popdata <- do.call("cbind.data.frame", popdatalist.new)
  
  names(popdata) <- c("rowid", "popid", "patchid", "individ", "year2", "firstseen", "lastseen", "obsage", 
                      "obslifespan", "xpos1", "ypos1", "sizea1", "sizeb1", "sizec1", "size1added", "repstra1", 
                      "repstrb1", "repstr1added", "feca1", "fecb1", "fec1added", "indcova1", "indcovb1",
                      "indcovc1", "censor1", "juvgiven1", "obsstatus1", "repstatus1", "fecstatus1", "matstatus1",
                      "alive1", "stage1", "stage1index", "xpos2", "ypos2", "sizea2", "sizeb2", "sizec2", "size2added",
                      "repstra2", "repstrb2", "repstr2added", "feca2", "fecb2", "fec2added", "indcova2", "indcovb2",
                      "indcovc2", "censor2", "juvgiven2", "obsstatus2", "repstatus2", "fecstatus2", "matstatus2", 
                      "alive2", "stage2", "stage2index", "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "size3added", 
                      "repstra3", "repstrb3", "repstr3added", "feca3", "fecb3", "fec3added", "indcova3", "indcovb3",
                      "indcovc3", "censor3", "juvgiven3", "obsstatus3", "repstatus3", "fecstatus3", "matstatus3", 
                      "alive3", "stage3", "stage3index")
  
  popdata <- subset(popdata, alive2 == 1)
  
  if (censor) {
    popdata <- subset(popdata, censor1 == censorkeep & censor2 == censorkeep)
    popdata <- subset(popdata, censor3 == censorkeep)
  }
  
  if (!is.na(spacing)) {
    popdata$density <- .density3(popdata, which(names(popdata) == "xpos2"), 
                                 which(names(popdata) == "ypos2"), 
                                 which(names(popdata) == "year2"), spacing)
  }
  
  if (reduce) {
    if (all(is.na(popdata$xpos1)) | length(unique(popdata$xpos1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "xpos1"))]
    } else if (all.equal(sort(unique(popdata$xpos1)), c(-1,0))) {popdata <- popdata[,-c(which(names(popdata) == "xpos1"))]}
    if (all(is.na(popdata$ypos1)) | length(unique(popdata$ypos1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "ypos1"))]
    } else if (all.equal(sort(unique(popdata$ypos1)), c(-1,0))) {popdata <- popdata[,-c(which(names(popdata) == "ypos1"))]}
    
    if (all(is.na(popdata$xpos2)) | length(unique(popdata$xpos2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "xpos2"))]
    } else if (all.equal(sort(unique(popdata$xpos2)), c(-1,0))) {popdata <- popdata[,-c(which(names(popdata) == "xpos2"))]}
    if (all(is.na(popdata$ypos2)) | length(unique(popdata$ypos2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "ypos2"))]
    } else if (all.equal(sort(unique(popdata$ypos2)), c(-1,0))) {popdata <- popdata[,-c(which(names(popdata) == "ypos2"))]}
    if (all(is.na(popdata$xpos3)) | length(unique(popdata$xpos3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "xpos3"))]
    } else if (all.equal(sort(unique(popdata$xpos3)), c(-1,0))) {popdata <- popdata[,-c(which(names(popdata) == "xpos3"))]}
    if (all(is.na(popdata$ypos3)) | length(unique(popdata$ypos3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "ypos3"))]
    } else if (all.equal(sort(unique(popdata$ypos3)), c(-1,0))) {popdata <- popdata[,-c(which(names(popdata) == "ypos3"))]}
    
    if (!censor) {
      popdata <- popdata[,-c(which(names(popdata) =="censor1"))]
      popdata <- popdata[,-c(which(names(popdata) =="censor2"))]
      popdata <- popdata[,-c(which(names(popdata) =="censor3"))]
    }
    
    if (all(is.na(popdata$sizea1)) | length(unique(popdata$sizea1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizea1"))]
    }
    if (all(is.na(popdata$sizea2)) | length(unique(popdata$sizea2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizea2"))]
    }
    if (all(is.na(popdata$sizea3)) | length(unique(popdata$sizea3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizea3"))]
    }
    
    if (all(is.na(popdata$sizeb1)) | length(unique(popdata$sizeb1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizeb1"))]
    }
    if (all(is.na(popdata$sizeb2)) | length(unique(popdata$sizeb2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizeb2"))]
    }
    if (all(is.na(popdata$sizeb3)) | length(unique(popdata$sizeb3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizeb3"))]
    }
    
    if (all(is.na(popdata$sizec1)) | length(unique(popdata$sizec1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizec1"))]
    }
    if (all(is.na(popdata$sizec2)) | length(unique(popdata$sizec2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizec2"))]
    }
    if (all(is.na(popdata$sizec3)) | length(unique(popdata$sizec3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="sizec3"))]
    }
    
    if (isTRUE(all.equal(popdata$size1added, popdata$sizea1))) {
      popdata <- popdata[,-c(which(names(popdata) == "size1added"))]
    } else if (all(is.na(popdata$size1added)) | length(unique(popdata$size1added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "size1added"))]
    }
    if (isTRUE(all.equal(popdata$size2added, popdata$sizea2))) {
      popdata <- popdata[,-c(which(names(popdata) == "size2added"))]
    } else if (all(is.na(popdata$size2added)) | length(unique(popdata$size2added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "size2added"))]
    }
    if (isTRUE(all.equal(popdata$size3added, popdata$sizea3))) {
      popdata <- popdata[,-c(which(names(popdata) == "size3added"))]
    } else if (all(is.na(popdata$size3added)) | length(unique(popdata$size3added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "size3added"))]
    }
    
    if (all(is.na(popdata$repstra1)) | length(unique(popdata$repstra1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "repstra1"))]
    }
    if (all(is.na(popdata$repstra2)) | length(unique(popdata$repstra2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "repstra2"))]
    }
    if (all(is.na(popdata$repstra3)) | length(unique(popdata$repstra3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "repstra3"))]
    }
    
    if (all(is.na(popdata$repstrb1)) | length(unique(popdata$repstrb1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "repstrb1"))]
    }
    if (all(is.na(popdata$repstrb2)) | length(unique(popdata$repstrb2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "repstrb2"))]
    }
    if (all(is.na(popdata$repstrb3)) | length(unique(popdata$repstrb3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) == "repstrb3"))]
    }
    
    if (isTRUE(all.equal(popdata$repstr1added, popdata$repstr1a))) {
      popdata <- popdata[,-c(which(names(popdata) =="repstr1added"))]
    } else if (all(is.na(popdata$repstr1added)) | length(unique(popdata$repstr1added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="repstr1added"))]
    }
    if (isTRUE(all.equal(popdata$repstr2added, popdata$repstr2a))) {
      popdata <- popdata[,-c(which(names(popdata) =="repstr2added"))]
    } else if (all(is.na(popdata$repstr2added)) | length(unique(popdata$repstr2added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="repstr2added"))]
    }
    if (isTRUE(all.equal(popdata$repstr3added, popdata$repstr3a))) {
      popdata <- popdata[,-c(which(names(popdata) =="repstr3added"))]
    } else if (all(is.na(popdata$repstr3added)) | length(unique(popdata$repstr3added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="repstr3added"))]
    }
    
    if (all(is.na(popdata$feca1)) | length(unique(popdata$feca1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="feca1"))]
    }
    if (all(is.na(popdata$feca2)) | length(unique(popdata$feca2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="feca2"))]
    }
    if (all(is.na(popdata$feca3)) | length(unique(popdata$feca3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="feca3"))]
    }
    
    if (all(is.na(popdata$fecb1)) | length(unique(popdata$fecb1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="fecb1"))]
    }
    if (all(is.na(popdata$fecb2)) | length(unique(popdata$fecb2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="fecb2"))]
    }
    if (all(is.na(popdata$fecb3)) | length(unique(popdata$fecb3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="fecb3"))]
    }
    
    if (isTRUE(all.equal(popdata$fec1added, popdata$feca1))) {
      popdata <- popdata[,-c(which(names(popdata) =="fec1added"))]
    } else if (all(is.na(popdata$fec1added)) | length(unique(popdata$fec1added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="fec1added"))]
    }
    if (isTRUE(all.equal(popdata$fec2added, popdata$feca2))) {
      popdata <- popdata[,-c(which(names(popdata) =="fec2added"))]
    } else if (all(is.na(popdata$fec2added)) | length(unique(popdata$fec2added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="fec2added"))]
    }
    if (isTRUE(all.equal(popdata$fec3added, popdata$feca3))) {
      popdata <- popdata[,-c(which(names(popdata) =="fec3added"))]
    } else if (all(is.na(popdata$fec3added)) | length(unique(popdata$fec3added)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="fec3added"))]
    }
    
    if (all(is.na(popdata$indcova1)) | length(unique(popdata$indcova1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcova1"))]
    }
    if (all(is.na(popdata$indcova2)) | length(unique(popdata$indcova2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcova2"))]
    }
    if (all(is.na(popdata$indcova3)) | length(unique(popdata$indcova3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcova3"))]
    }
    if (all(is.na(popdata$indcovb1)) | length(unique(popdata$indcovb1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcovb1"))]
    }
    if (all(is.na(popdata$indcovb2)) | length(unique(popdata$indcovb2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcovb2"))]
    }
    if (all(is.na(popdata$indcovb3)) | length(unique(popdata$indcovb3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcovb3"))]
    }
    if (all(is.na(popdata$indcovc1)) | length(unique(popdata$indcovc1)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcovc1"))]
    }
    if (all(is.na(popdata$indcovc2)) | length(unique(popdata$indcovc2)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcovc2"))]
    }
    if (all(is.na(popdata$indcovc3)) | length(unique(popdata$indcovc3)) == 1) {
      popdata <- popdata[,-c(which(names(popdata) =="indcovc3"))]
    }
    
    if (all(popdata$obsstatus1 == popdata$obsstatus1[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="obsstatus1"))]
    }
    if (all(popdata$obsstatus2 == popdata$obsstatus2[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="obsstatus2"))]
    }
    if (all(popdata$obsstatus3 == popdata$obsstatus3[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="obsstatus3"))]
    }
    
    if (all(popdata$repstatus1 == popdata$repstatus1[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="repstatus1"))]
    }
    if (all(popdata$repstatus2 == popdata$repstatus2[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="repstatus2"))]
    }
    if (all(popdata$repstatus3 == popdata$repstatus3[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="repstatus3"))]
    }
    
    if (all(popdata$fecstatus1 == popdata$fecstatus1[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="fecstatus1"))]
    }
    if (all(popdata$fecstatus2 == popdata$fecstatus2[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="fecstatus2"))]
    }
    if (all(popdata$fecstatus3 == popdata$fecstatus3[1])) {
      popdata <- popdata[,-c(which(names(popdata) =="fecstatus3"))]
    }
  }
  
  class(popdata) <- append(class(popdata), "hfvdata")
  
  return(popdata)
}

#' Estimate Density on Basis of Cartesian Coordinates
#' 
#' \code{.density3()} estimates density on the basis of Cartesian coordinates and
#' spacing information supplied as input. It is used internally by 
#' \code{\link{historicalize3}()} and \code{\link{verticalize3}()}.
#' 
#' @param data Demographic dataset in historical vertical format.
#' @param xcol Number of column in \code{data} corresponding to x position.
#' @param ycol Number of column in \code{data} corresponding to y position.
#' @param yearcol Number of column in \code{data} corresponding to time.
#' @param spacing Resolution of density estimation, as a scalar numeric.
#' 
#' @return This function returns the original data frame supplied as \code{data} but
#' with a new variable added to the end of the data frame showing the local density
#' of the individual.
#' 
#' @keywords internal
#' @noRd
.density3 <- function(data, xcol, ycol, yearcol, spacing) {
  data$Xgen <- floor(data[,xcol] / spacing)
  data$Ygen <- floor(data[,ycol] / spacing)
  
  data$grid.id <- paste(data$Xgen, data$Ygen, sep = ",")
  
  grid.densities <- xtabs(paste("~ grid.id + ", yearcol), data = data)
  
  data$density.est <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
    grid.densities[as.character(data$grid.id[X]), as.character(data[X, yearcol])]
  })
  
  data$density.est[which(is.na(data$Xgen))] <- NA
  
  return(data)
}
