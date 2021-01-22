plot.ictest <- function (x, which="all", ...)
{
 which <- match.arg(which, c("all", "k"))
 if (which == "all") {S <- x$S} else {S <- x$S[ ,0:x$k, drop=FALSE]}
 if(ncol(S) <= 2) {plot(S, ...)} else {
        if (any(class(S) %in% c("mts", "xts", "zoo"))) {plot(S, ...)} else {pairs(S, ...)}
        }
}
