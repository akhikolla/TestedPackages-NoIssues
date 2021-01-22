check <- modules::module({

  samplingVar <- function(samplingVar) {

    anyNa <- any(is.na(samplingVar))

    if (anyNa) stop("Missing values in sampling variances are not allowed.")

    anyZero <- any(samplingVar <= 0)

    if (anyZero) {
      warning("There are zero sampling variances. They are replaced by 1e-05.")
      samplingVar <- replace(samplingVar, samplingVar == 0, 1e-05)
    }

    samplingVar
    
  }
  
})
