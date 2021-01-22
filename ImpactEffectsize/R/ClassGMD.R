#Pools two values of Gini's mean difference calculated for separate groups.
#' @importFrom stats var
# gini.md<- function(x) {
#   j <-order(x)
#   n <-length(x)
#   return(4*sum((1:length(x))*x[j])/(n*(n-1))-2*mean(x)*(n+1)/(n-1))
# }

ClassGMD <- function(Data, Cls) {
  maxPoints = 10000
  if (length(Data) > maxPoints) {
    df <- data.frame(cbind(Data, Cls))
    table(df$Cls)
    dfsplit <- split(df, list(df$Cls))
    set.seed(42)
    samples <-
      lapply(dfsplit, function(x)
        x[sample(1:nrow(x), maxPoints, FALSE),])
    out <- do.call(rbind, samples)
    table(out$Cls)
    Data <- as.vector(out$Data)
    Cls <- as.vector(out$Cls)
  }
  if (var(Data) == 0)
    GMDn <- 1e-7
  else if ((var(Data[Cls == sort(unique(Cls))[1]]) == 0 |
            var(Data[Cls == sort(unique(Cls))[2]]) == 0) &
           var(Data) > 0)
    GMDn <- c_gmd(Data)
  else {
    GMD1 <- c_gmd(Data[Cls == sort(unique(Cls))[1]])
    GMD2 <- c_gmd(Data[Cls == sort(unique(Cls))[2]])
    GMDn <- sqrt((GMD1 ^ 2 + GMD2 ^ 2) / 2)
  }
  return(GMDn)
}
