cautiousLearning <- function(chart, x, mu0, s0, plot=TRUE) {
    y <- ts(applyChart(chart, x, mu0, s0))
    if (chart$chart == "X") {
        nm <- "X"
    } else if (chart$chart == "EWMA") {
        nm <- "EWMA"
    }  else {
        nm <- c("C-", "C+")
    }
    colnames(y) <- c(nm, "CL", "LCL", "UCL", "MUhat", "Shat", "Lhat")
    if (plot) {
        nstat <- NCOL(y)-6
        r <- range(y[,1:(nstat+3)])
        d <- (r[2]-r[1])/10
        plot.ts(y[,1:nstat], plot.type="single", type="b",
                ylab=chart$chart, ylim=r+c(-d,+d))
        lines(y[,nstat+1], type="s", lty="solid",  col="red")
        lines(y[,nstat+2], type="s", lty="dashed", col="red")
        lines(y[,nstat+3], type="s", lty="dashed", col="red")
        invisible(y)
    } else {
        y
    }
}
