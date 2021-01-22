print.ypinterim <- function(x, ...) {

    digit <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")
    x_teststat <- sprintf(digit, x$teststat)
    x_pvalue <- fun_less(sprintf(digit, x$pvalue))
    x_pvalue_bnd <- fun_less(sprintf(digit, x$pvalue.bnd))
    x_critvalue <- sprintf(digit, x$all.critvalue)
    x_usercritvalue <- sprintf(digit, x$user.critvalue)
    x_nlooks <- x$n.looks
    x_bound <- x$bound
    x_repnum <- x$repnum
    x_names <- paste("Look", 1:x_nlooks)
    x_spendfun <- fun_less(sprintf(digit, x$dspendfun))
    cx_spendfun <- fun_less(sprintf(digit, x$spendfun))

    cat("\n===========================================================\n")
    cat("       Interim Testing For Survival Outcomes")
    cat("\n===========================================================\n")

    name_testp <- rep("       TestStat (nominal p-value):", length(x_names))
    name_testp[1] <- "       TestStat (nominal p-value*):"
    name_crtp <- rep("       Boundary (nominal p-value):", length(x_names))
    
    for (ilook in 1:length(x_names)) {
        cat("    <", x_names[ilook], ">", "\n\n")
        cat("       Alpha allocated (cumulative):", paste(x_spendfun[ilook]," (", cx_spendfun[ilook], ")", sep=""), "\n")
        cat(name_testp[ilook], x_teststat[ilook], paste("(", x_pvalue[ilook],
            ")", sep = ""), "\n")
        cat(name_crtp[ilook], x_critvalue[ilook], paste("(", x_pvalue_bnd[ilook],
                                                        ")", sep = ""), "\n")
        if (ilook != length(x_names))
            cat("\n-------------------------------------------------------\n")
    }
    cat("\n===========================================================\n")
    cat("* The probability of the test statistic exceeding the critical value\n")
    cat(" at that look, regardless of the test behavior at other looks.\n")
}

