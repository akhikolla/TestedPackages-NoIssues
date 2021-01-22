screeplot.ictest <- function (x, type = "barplot", main = deparse(substitute(x)),
    ylab = "criterion", xlab = "component", ...)
{
    type <- match.arg(type, c("barplot", "lines"))

    if (type == "barplot") {
        barplot(x$D, ylab = ylab, xlab = xlab, main = main, ...)
    } else {
        index <- 1:length(x$D)
        plot(index, x$D, type = "b", ylab = ylab,
            xlab = xlab, axes = FALSE, main = main, ...)
        axis(2)
        axis(1, at = index, labels = TRUE)
    }
}
