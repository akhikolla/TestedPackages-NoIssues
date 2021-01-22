\name{incomeUK}
\alias{incomeUK}
\usage{data(incomeUK)}
\title{Incomes of United Kingdom citizens}
\description{
The \code{incomeUK} numeric vector has 7,201 incomes of United Kingdom citizens for the year 1975. The data have been divided by average income.
}
\format{
Each entry of the numeric vector is the ratio of the income of a British citizen and the average income of the sample.
}
\source{The Economic and Social Research Council Data Archive at the University of Essex, United Kingdom, Family Expenditure Survey, Annual Tapes, 1968-1983, Department of Employment, Statistics Division, Her Majesty's Stationery Office, London.
}
\examples{
library(densEstBayes)
data(incomeUK) 
hist(incomeUK,breaks = 100,col = "gold")

# Obtain and plot ordinary density estimate:

dest <- densEstBayes(incomeUK,method = "SMFVB",
                     control = densEstBayes.control(range.x = c(0,4)))
plot(dest,xlab = "income (multiple of average income)")
rug(incomeUK,col = "dodgerblue",quiet = TRUE)

# Now obtain and plot improved density estimate using log transformation:

destlogScale <- densEstBayes(log(incomeUK),method = "SMFVB",
                          control = densEstBayes.control(range.x = c(-3,1.2)))
plotVecs <- plot(destlogScale,plotIt = FALSE)
xLogScaleg <- plotVecs$xg
densLowLogScaleg <- plotVecs$densLowg
densEstLogScaleg <- plotVecs$densEstg
densUppLogScaleg <- plotVecs$densUppg
xg <- exp(xLogScaleg)
densLowg <- densLowLogScaleg/xg
densEstg <- densEstLogScaleg/xg
densUppg <- densUppLogScaleg/xg
plot(0,xlim = range(xg),ylim = range(c(densLowg,densUppg)),type = "n",bty = "l",
     xlab = "income (multiple of average income)",ylab = "density")
polygon(c(xg,rev(xg)),c(densLowg,rev(densUppg)),
        col="palegreen",border = FALSE)
lines(xg,densEstg,col = "darkgreen",lwd = 2)
abline(0,0,col = "slateblue")
rug(incomeUK,col = "dodgerblue",quiet = TRUE)
}


