# Take data set in my usual form and return a data set under LVCF.
# It returns it in the form of a data.frame called new.data
# @title FUNCTION_TITLE
# @description FUNCTION_DESCRIPTION
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
# @param obs.tm Vector of observed main event time or censoring time
# @param delta Vector of censoring indicators. \code{1} for event \code{0} for censored
# @param Z Additional variables for the main model other than the binary covaraite, Default: NULL
# @return OUTPUT_DESCRIPTION
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @rdname LVCFdata
# @export 
LVCFdata <- function(w, w.res, obs.tm, delta, Z = NULL)
{
  n.sample <- length(obs.tm)
  if (is.null(Z))
  {
  new.data <- matrix(nrow = 2*n.sample, ncol =  5) # 5:  1 for ID,  2 for (start,stop) time, one for X and one for delta
  colnames(new.data) <- c("ID", "start.time", "stop.time", "delta", "X")
  } else {
  new.data <- matrix(nrow = 2*n.sample, ncol = ncol(Z) + 5) # 5:  1 for ID,  2 for (start,stop) time, one for X and one for delta
  colnames(new.data) <- c("ID", "start.time", "stop.time", "delta", "X", paste0("Z", 1:ncol(Z)))
  }
  first.one <- apply(w.res, 1, function(x) Position(function(y) y==1, x))  # Find for each observation the first time w.res==1
  k <- 1
  # For each observation, set one row if X(t)==0 for all observed t, and set two rows (one with X==0 and one with X==1) if a change was observed
  for (j in 1:n.sample)
  {
    if (is.na(first.one[j]))
    {
      if (is.null(Z))
      {
        new.data[k, ] <- c(j, 0, obs.tm[j], delta[j], 0) 
      } else { 
    new.data[k, ] <- c(j, 0, obs.tm[j], delta[j], 0, Z[j,]) 
      }
    k <- k + 1
    } else 
      {
        change.point <- w[j, first.one[j]] # when X(t)==1 was first observed
        if (change.point > obs.tm[j]) {
          if (is.null(Z))
          { new.data[k, ] <- c(j, 0, obs.tm[j], delta[j], 0) 
          } else {
          new.data[k, ] <- c(j, 0, obs.tm[j], delta[j], 0, Z[j,]) 
          }
          k <- k + 1
          } else {
            if (is.null(Z))
            {
              new.data[k, ] <- c(j, 0, change.point, 0, 0)
              new.data[k + 1, ] <- c(j, change.point, obs.tm[j], delta[j], 1)
            } else {
        new.data[k, ] <- c(j, 0, change.point, 0, 0, Z[j, ])
        new.data[k + 1, ] <- c(j, change.point, obs.tm[j], delta[j], 1, Z[j, ])
            }
        k <- k + 2
      }}
  }
  new.data <- as.data.frame(new.data[1:(k - 1), ]) # k + 1 was the last row, but then we added +2 to k
  return(as.data.frame(new.data)) # return a data.frame
}
