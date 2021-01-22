covSIR <- function(X, y, h = 10, ...) {
  n <- nrow(X)
  slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/h), ...)),
                          include.lowest = TRUE ,labels = FALSE))
  slicemean <- do.call(rbind, by(X, slices, colMeans))
  nh <- as.numeric(table(slices))
  cov.wt(slicemean, wt=nh/n, method="ML")$cov
  #cov(slicemean)
}
