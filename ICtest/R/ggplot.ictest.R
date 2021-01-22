ggplot.ictest <- function (data, mapping=aes(), mapvar=NULL, which="all", ..., environment=parent.frame())
{
 x <- data
 which <- match.arg(which, c("all", "k"))
 if (which=="all") {S <- x$S} else {S <- x$S[ ,0:x$k, drop=FALSE]}
 
 if(!is.null(mapvar)){
   dummy <- data.frame(S, mapvar)
 }else{
   dummy <- data.frame(S)
 }
 original_columns <- colnames(S)
 if (any(class(S) %in% c("mts", "xts", "zoo"))) {
   p <- ncol(S)
   n <- nrow(S)
   S_long <- data.frame(var = c(as.matrix(S)), lat = factor(rep(1:p, each = n), labels = sapply(1:p, function(x) paste0("Series ", x))), time = rep(1:n, times = p))
   ggplot(S_long, aes(x = time, y = var)) +
     geom_line() +
     facet_wrap(~ lat, scales = "free_y", ncol = ceiling((p + 1)/5)) +
     labs(x = "Time", y = "")
 } else {
   ggpairs(dummy, mapping=mapping, columns = original_columns, ...)
 }
}
