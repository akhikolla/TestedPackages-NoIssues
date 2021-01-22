# function for a ggplot ladleplot - user can choose the criterion
ggladleplot <- function(x, crit = "gn", type="l", ylab = crit, xlab = "component", main = deparse(substitute(x)), ...)
{
  crit <- match.arg(crit, c("gn", "fn", "phin", "lambda"))
  index <- 0:(length(x[[crit]])-1)
  plot_data <- data.frame(crit = x[[crit]], index = index)
  ggplot(plot_data, aes(x = index, y = crit)) +
    geom_line() +
    geom_point() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(main)
}

