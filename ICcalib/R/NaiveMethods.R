#### Functions for the naive methods that do not perform calibration###
# A new version is LVCFdata and MidIdata

# All three functions take three arguments:
# w: follow-up times, where infinite time implies no more follow up
# w.res: reult in each follow-up time: zero or one.
# point: the point in which a value for X is needed (i.e., current risk.set)
################################################################################################################
################### CarryForward ###############################################################################
CarryForward <- function(w, w.res, point) {
  n.sample <- nrow(w)
  interval.w <-   FindIntervalCPP(point = point, w =w)
  #  interval.w[is.na(interval.w)] <- ncol(w) + 1 
  x.carry <- vector(length = n.sample)
  for (j in 1:n.sample)
  {
    if (interval.w[j]==1)
    {
      x.carry[j] <- 0
    } else  if (any(w.res[j,1:(interval.w[j]-1)]==1))
    {
      x.carry[j] <- 1
    } else
    {
      x.carry[j] <- 0
    }} 
  return(x.carry)
  return(interval.w)
}
################################################################################################################
################### CarryBack ##################################################################################
CarryBack <- function(w, w.res, point) {
  n.sample <- nrow(w)
  n.intervals <- ncol(w)
  interval.w <-   FindIntervalCPP(point = point, w =w)
  #  interval.w[is.na(interval.w)] <- ncol(w) + 1 
  x.carry <- vector(length = n.sample)
  for (j in 1:n.sample)
  {
    if (interval.w[j]==1)
    {
      if(w.res[j,1]==1)
      {       x.carry[j] <- 1} else {x.carry[j] <- 0}
    } else if (interval.w[j]==n.intervals+1)
    {
      if(any(w.res[j,1:n.intervals]==1))
      {       x.carry[j] <- 1} else {x.carry[j] <- 0}    
    }    
    else  if (any(w.res[j,1:interval.w[j]]==1))
    {
      x.carry[j] <- 1
    } else
    {
      x.carry[j] <- 0
    }} 
  return(x.carry)
}
################################################################################################################
################### Midpoint ###################################################################################
MidPoint <- function(w, w.res, point) 
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  colnames(lr.for.fit) <- c("left","right")
  tm.v <- vector(length=nrow(w)) #time of chhange
  tm.v[lr.for.fit$right==Inf] <- Inf
  tm.v[lr.for.fit$right<Inf] <- (lr.for.fit$right[lr.for.fit$right<Inf]+lr.for.fit$left[lr.for.fit$right<Inf])/2
  x <- ifelse(tm.v<point,1,0)
  return(x)
}
