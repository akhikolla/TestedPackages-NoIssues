summary.ypinterim <- function(object, ...) {
    x <- object
    if (class(x) != "ypinterim")
        stop("Please use the object from the 'ypinterim' function")

    digit <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")
    x_teststat <- sprintf(digit, x$teststat)
    x_pvalue <- fun_less(sprintf(digit, x$pvalue))
    x_critvalue <- sprintf(digit, x$all.critvalue)
    x_pvalue_bound <- fun_less(sprintf(digit, x$pvalue.bnd))
    x_bound <- x$bound
    x_repnum <- x$repnum
    x_length <- length(x_teststat)
    x_names <- paste("Look", 1:x_length)
    x_spendfun <- fun_less(sprintf(digit, x$dspendfun))
    cx_spendfun <- fun_less(sprintf(digit, x$spendfun))

    test_and_p <- paste(x_teststat, " (", x_pvalue, ")", sep = "")
    bnd_and_p <- paste(x_critvalue, " (", x_pvalue_bound, ")", sep = "")
    sresult <- cbind(x_spendfun, cx_spendfun, test_and_p, bnd_and_p)
    
    colnames(sresult) <- c("alpha.spent","cum.alpha", "test.stat (nominal p)", "boundary (nominal p)")
    rownames(sresult) <- paste("Look", "(", 1:x_length, ")", sep = "")

    as.data.frame(sresult)
}
