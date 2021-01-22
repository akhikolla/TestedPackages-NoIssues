## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(tvR)

## ----example1------------------------------------------------------------
set.seed(1)
x = rep(sample(1:5,10,replace=TRUE), each=50) ## main signal
xnoised = x + rnorm(length(x), sd=0.25)       ## add noise

## ----compare1------------------------------------------------------------
## apply denoising process
xproc1 = denoise1(xnoised, method = "TVL2.IC")
xproc2 = denoise1(xnoised, method = "TVL2.MM")

## ----compare1vis, echo=FALSE, fig.align='center', fig.fullwidth=TRUE, fig.width=7, fig.height=5----
## plot noisy and denoised signals
plot(xnoised, pch=19, cex=0.1, main="compare two algorithms", xlab="time domain", ylab="signal value")
lines(xproc1, col="blue", lwd=2)
lines(xproc2, col="red", lwd=2)
legend("topright",legend=c("Noisy","TVL2.IC","TVL2.MM"),
col=c("black","blue","red"),#' lty = c("solid", "solid", "solid"),
lwd = c(0, 2, 2), pch = c(19, NA, NA),
pt.cex = c(1, NA, NA), inset = 0.05)

## ----compare2------------------------------------------------------------
compare = list()
for (i in 1:4){
  compare[[i]] = denoise1(xnoised, lambda = 10^(i-4), method="TVL2.IC")
}

## ----compare2vis, fig.show='hold', echo=FALSE----------------------------
for (i in 1:4){
  pm = paste("lambda=1e",i-4,sep="")
  plot(xnoised, pch=19, cex=0.2, main=pm, xlab="time domain", ylab="signal value")
  lines(compare[[i]], col=as.integer(i+2), lwd=1.2)
}

## ----lena----------------------------------------------------------------
data(lena128)
xnoised <- lena128 + array(rnorm(128*128, sd=10), c(128,128))

## ----compare3------------------------------------------------------------
## apply denoising process
xproc1 <- denoise2(xnoised, lambda=10, method="TVL1.PrimalDual")
xproc2 <- denoise2(xnoised, lambda=10, method="TVL2.FiniteDifference")
xproc3 <- denoise2(xnoised, lambda=10, method="TVL2.PrimalDual")

## ----compare3vis, echo=FALSE, fig.show='hold'----------------------------
gcol = gray(0:128/128)
image(xnoised, main="Noised", col=gcol)
image(xproc1, main="L1-PrimalDual", col=gcol)
image(xproc2, main="L2-FiniteDifference", col=gcol)
image(xproc3, main="L2-PrimalDual", col=gcol)

