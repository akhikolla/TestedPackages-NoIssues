
getCI = function(t, ts, a){
  times = length(ts)
    #Normal CI
    z975 = qnorm(1-a/2, 0,1)
    NorCI = c(t-z975*sqrt(var(ts)), t+z975*sqrt(var(ts)))

    #Pivotal CI
    PivCI = c(t-sort(ts)[(1-a/2)*times], t-sort(ts)[(a/2)*times])
    if(t < PivCI[1] || t > PivCI[2])
      warning("Pivotal CI doesn't include test statistic, result might be misleading.")

    #Percentile CI
    PerCI = c(sort(ts)[(a/2)*times], sort(ts)[(1-a/2)*times])

    CI = list(NorCI, PivCI, PerCI)
    names(CI) = c("Normal CI", "Pivotal CI", "Percentile CI")
    return (CI)
}
