### PREDICTION ####
methods::setMethod(
  f="predict",
  signature = "ResultClassifOrdinal",
  definition = function(object, x=matrix(0,nrow=1,ncol=1)) {
    res <- predictions(object, x)
    return(res)
  }
)

### PLOTTING ###
methods::setGeneric("plot",function(object){standardGeneric("plot")}) 

methods::setMethod(
  f="plot",
  signature = c("ResultClassifOrdinal"),
  definition = function(object) {
    bosplot(object)
  }
)

methods::setMethod(
  f="plot",
  signature = c("ResultCoclustOrdinal"),
  definition = function(object) {
    bosplot(object)
  }
)

methods::setMethod(
  f="plot",
  signature = c("ResultClustOrdinal"),
  definition = function(object) {
    bosplot(object)
  }
)

### summaryING ###
methods::setGeneric("summary",function(object){standardGeneric("summary")}) 

methods::setMethod(
  f="summary",
  signature = c("ResultClassifOrdinal"),
  definition = function(object) {
    bossummary(object)
  }
)

methods::setMethod(
  f="summary",
  signature = c("ResultCoclustOrdinal"),
  definition = function(object) {
    bossummary(object)
  }
)

methods::setMethod(
  f="summary",
  signature = c("ResultClustOrdinal"),
  definition = function(object) {
    bossummary(object)
  }
)


