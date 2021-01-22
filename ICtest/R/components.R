components <- function(x, ...) UseMethod("components")

components.ictest <- function (x, which="all", ...)
{
 which <- match.arg(which, c("all", "k"))
 if (which=="all") {S <- x$S} else {S <- x$S[ ,0:x$k, drop=FALSE]}
 S
}
