# Estimates the Pareto density function (PDE)
#  Author: MT 2019 (modified)
#' @importFrom stats, runif, sd, quantile, quantile
OptimalNoBinsIE <- function(Data) {
  if (is.matrix(Data))
    nData <- colSums(!is.nan(Data))
  if (is.vector(Data))
    nData <- sum(!is.nan(Data))
  
  if (nData < 1) {
    optNrOfBins <- 0
  } else{
    sigma <- sd(Data, na.rm = TRUE)
    ifelse(nData < 5000, p <- c_quantile(na.omit(Data), c(0.25, 0.75)), p <- quantile(Data, c(0.25, 0.75), type = 8, na.rm = TRUE)    )
    interquartilRange <- p[2] - p[1]
    sigmaSir <- min(sigma, interquartilRange / 1.349)
    optBinWidth <- 3.49 * sigmaSir / (nData) ^ (1 / 3)
    if (optBinWidth > 0) {
      optNrOfBins <-  max(ceiling((max(Data, na.rm = TRUE) - min(Data, na.rm = TRUE)) / optBinWidth), 10)
    } else{
      optNrOfBins <- 10
    }
  }
  return (optNrOfBins)
}


ParetoRadiusIE <- function(Data,  maximumNrSamples = 10000) {
  requireNamespace('parallelDist')
  ntemp <- sum(is.nan(Data))
  if (ntemp > 0) {
    warning('Data has NaN values, Pareto Radius may not be calculated.')
  }
  ntemp2 <- sum(is.na(Data))
  if (ntemp2 > ntemp)
    warning('Data has NA values, Pareto Radius may not be calculated.')
  ntemp3 <- sum(is.infinite(Data))
  if (ntemp3 > 0)
    warning('Data has infinite valuues, Pareto Radius may not be calculated.')
  
  nData <- length(Data)
  
  if (maximumNrSamples >= nData) {
    # no sampling necessary
    sampleData <- Data
  } else{
    #  sample with uniform distribution MaximumNrSamples
    sampleInd <- ceiling(runif(maximumNrSamples, min = 0, max = nData)) # floor(nData*c(runif(maximumNrSamples))+1)
    sampleData <- Data[sampleInd]
  }
  # calculate distances
  distvec <-
    as.vector(parallelDist::parallelDist(as.matrix(sampleData), method = 'euclidean', upper = FALSE, diag = FALSE))
  # selection of ParetoRadius
  ifelse(nData < 5000, paretoRadius <- c_quantile(na.omit(distvec), 18 / 100),
    paretoRadius <- quantile(distvec, 18 / 100, type = 8, na.rm = TRUE))
  
  
  if (paretoRadius == 0)
  {
    ifelse(nData < 5000, pzt <- c_quantile(na.omit(distvec), probs = c(1:100) / 100),
      pzt <- quantile(distvec, probs = c(1:100) / 100, type = 8, na.rm = TRUE)
    )
    paretoRadius <- min(pzt[pzt > 0], na.rm = TRUE) # take the smallest nonzero
  }
  if (is.nan(paretoRadius))
    stop('Pareto Radius could not be calculated. (NaN value)')
  if (is.na(paretoRadius))
    stop('Pareto Radius could not be calculated. (NA value)')
  if (!is.finite(paretoRadius))
    stop('Pareto Radius could not be calculated. (infinite value)')
  #    plot of distance distribution
  if (nData > 1024)
    paretoRadius = paretoRadius * 4 / (nData ^ 0.2)
  return(paretoRadius)
}


ParetoDensityEstimationIE <- function(Data, paretoRadius, kernels = NULL, MinAnzKernels = 100) {
    requireNamespace('caTools') #fuer trapz
    if (!is.vector(Data)) {
      Data <- as.vector(Data)
      warning('Beware: ParetoDensityEstimationIE: Data set not univariate !')
    }
    if (!is.numeric(Data)) {
      Data <- as.numeric(Data)
      warning('Beware: ParetoDensityEstimationIE: Data set not numeric !')
    }
    if (length(Data) != sum(is.finite(Data))) {
      message('Not all values are finite. Please check of infinite or missing values.')
    }
    Data <- Data[is.finite(Data)]
    values <- unique(Data)
    
    if (length(values) > 2 & length(values) < 5) {
      warning('Less than 5 unqiue values for density estimation. Function may not work')
    }
    if (length(values) < 3) {
      warning(
        '1 or 2 unique values for density estimation. Dirac Delta distribution(s) is(are) assumed. Input of "kernels", "paretoRadius" and "MinAnzKernels" or ignored!'
      )
      if (values[1] != 0)
        kernels <- seq(from = values[1] * 0.9,
                       to = values[1] * 1.1,
                       by = values[1] * 0.0001)
      else
        kernels <- seq(from = values[1] - 0.1,
                       to = values[1] + 0.1,
                       by = 0.0001)
      
      paretoDensity <- rep(0, length(kernels))
      paretoDensity[kernels == values[1]] = 1
      
      if (length(values) == 2) {
        if (values[2] != 0)
          kernels2 <- seq(from = values[2] * 0.9,
                          to = values[2] * 1.1,
                          by = values[2] * 0.0001)
        else
          kernels2 <- seq(from = values[2] - 0.1,
                          to = values[2] + 0.1,
                          by = 0.0001)
        
        
        paretoDensity2 = rep(0, length(kernels2))
        paretoDensity2[kernels2 == values[2]] = 1
        
        paretoDensity <- c(paretoDensity, paretoDensity2)
        kernels = c(kernels, kernels2)
      }
      return (list(
        kernels = kernels,
        paretoDensity = paretoDensity,
        paretoRadius = 0
      ))
    }
    
    if (length(Data) < 10) {
      warning('Less than 10 datapoints given, ParetoRadiusIE potientially cannot be calcualted.')
    }
    if (missing(paretoRadius)) {
      paretoRadius <- ParetoRadiusIE(Data)
    } else if (is.null(paretoRadius)) {
      paretoRadius <- ParetoRadiusIE(Data)
    } else if (is.na(paretoRadius)) {
      paretoRadius <- ParetoRadiusIE(Data)
    } else if (paretoRadius == 0 || length(paretoRadius) == 0) {
      paretoRadius <- ParetoRadiusIE(Data)
    }
    
    minData <- min(Data, na.rm = TRUE)
    maxData <- max(Data, na.rm = TRUE)
    if (length(kernels) <= 1) {
      if (length(kernels) == 0 || (length(kernels) == 1 & kernels == 0)) {
        #MT: Korrektur: statt kernels==0 und im Input Kernels=0
        nBins <- OptimalNoBinsIE(Data)
        #MT: MinAnzKernels fehlte
        nBins <- max(MinAnzKernels , nBins)
        # mindestzahl von Kernels sicherstellen
        if (nBins > 100) {
          if (nBins > 1E4) {
            #MT: Fehlerabdfang bei zu vielen Bins
            nBins = 1E4
            warning('Too many bins estimated, try to transform or sample the data')
          } else{
            nBins = nBins * 3 + 1
          }
        }
        breaks <- pretty(c(minData, maxData), n = nBins, min.n = 1)
        nB <- length(breaks)
        mids <- 0.5 * (breaks[-1L] + breaks[-nB])
        kernels <- mids
      }
    }
    nKernels <- length(kernels)
    lowBInd <-  (Data < (minData + paretoRadius))
    lowR <- as.matrix(2 * minData - Data[lowBInd], ncol = 1)
    upBInd <-  (Data > (maxData - paretoRadius))
    upR <- as.matrix(2 * maxData - Data[upBInd], ncol = 1)
    DataPlus <- as.matrix(c(Data, lowR, upR), 1)
    
    paretoDensity <- rep(0, nKernels)
    for (i in 1:nKernels) {
      lb <- kernels[i] - paretoRadius
      ub <- kernels[i] + paretoRadius
      isInParetoSphere <- (DataPlus >= lb) & (DataPlus <= ub)
      paretoDensity[i] <- sum(isInParetoSphere)
    }
    
    area <- caTools::trapz(kernels, paretoDensity)
    if (area < 0.0000000001 || is.na(area)) {
      paretoDensity <- rep(0, nKernels)
    } else{
      paretoDensity <- paretoDensity / area
    }
    return (list(
      kernels = kernels,
      paretoDensity = paretoDensity,
      paretoRadius = paretoRadius
    ))
  }
