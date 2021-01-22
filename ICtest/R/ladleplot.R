# function for a ladleplot - user can choose the criterion
# ALSO a ggplot version? Like ggladleplot?
ladleplot <- function(x, crit = "gn", type="l", ylab = crit, xlab = "component", main = deparse(substitute(x)), ...)
{
  crit <- match.arg(crit, c("gn", "fn", "phin", "lambda"))
  index <- 0:(length(x[[crit]])-1)
  plot(index, x[[crit]], type = type, xlab = xlab, ylab=ylab, main = main, axes = FALSE, ...)
  axis(2)
  axis(1, at = index, labels = TRUE)
}
