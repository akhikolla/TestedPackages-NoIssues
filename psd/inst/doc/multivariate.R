## ----eval=TRUE, echo=FALSE----------------------------------------------------
rm(list=ls())
rng_et <- c(0.7, 2.2)
rng_ba <- c(0.0, 10)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
library(psd)
data(wipp30)
str(wipp30)
dat <- wipp30[, c(3,4,2)] # output as last column, input as first two
dat <- window(ts(dat), 100, 2400)
freqs <- 1:(nrow(dat)/2) * (24) / nrow(dat) # frequency in cpd

## ----eval=TRUE, echo=TRUE, label="calcmv"-------------------------------------
mv <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE, plot = FALSE)

## ----eval=TRUE, echo=TRUE, label="names"--------------------------------------
names(mv)

## ----eval=TRUE, echo=TRUE, label="plotnames", fig.height = 7------------------
pspectrum(dat, riedsid_column= 0L, verbose = FALSE, plot = TRUE)

## ----eval=TRUE, echo=TRUE, label="Tapers"-------------------------------------

mn <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE)   # minimum number 
mx <- pspectrum(dat, riedsid_column=-1L, verbose = FALSE)   # maximum number
c1 <- pspectrum(dat, riedsid_column= 1L, verbose = FALSE)   # column 1
c2 <- pspectrum(dat, riedsid_column= 2L, verbose = FALSE)   # column 2
c3 <- pspectrum(dat, riedsid_column= 3L, verbose = FALSE)   # column 3


## ----eval=TRUE, echo=FALSE, label="Taper plot", fig.height = 4.8--------------

ylim <- c(0, max(mx$taper))
layout(matrix(1:6, ncol = 1, nrow = 6), heights = c(rep(0.19, 5), 0.05))
par(mar = c(0.5, 4, 0, 0), mgp = c(2, 0.5, 0))
plot(mn$taper, ylim = ylim, xaxt="n")
legend('topleft', 'minimum')
plot(mx$taper, ylim = ylim, xaxt="n")
legend('topleft', 'maximum')
plot(c1$taper, ylim = ylim, xaxt="n")
legend('topleft', 'Water level column')
plot(c2$taper, ylim = ylim, xaxt="n")
legend('topleft', 'Barometric column')
plot(c3$taper, ylim = ylim)
legend('topleft', 'Earth tide column')
mtext('taper index', 1, line = 1.2, cex = 0.8)


## ----eval=TRUE, echo=TRUE, label="taper specification", fig.height = 3--------
par(mar = c(3, 4, 1, 0), mgp = c(2, 0.5, 0))
one_val <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE,
                 niter = 0, ntap.init = 11, plot = TRUE)
plot(one_val$taper)

## ----eval=TRUE, echo=TRUE, label="taper specification2", fig.height = 3-------
par(mar = c(3, 4, 1, 0), mgp = c(2, 0.5, 0))
vec_val <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE,
                 niter = 0, ntap = as.tapers(1:(nrow(dat) %/% 2) %/% 7 + 7),
                 plot = TRUE)
plot(vec_val$taper)

## ----eval=TRUE, echo=TRUE, label="pspec"--------------------------------------

spec <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE)$pspec
dim(spec)


## ----eval=TRUE, echo=TRUE, label="coh"----------------------------------------
coh <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE,
                 niter = 0, ntap.init = 11)$coh
coh_base <- spec.pgram(dat[, c(3,1,2)], spans = c(11), fast = FALSE, plot = FALSE)$coh

## ----eval=TRUE, echo=FALSE, label="coh_ba_plot", fig.height=3.3---------------
par(mar = c(3, 3, 2, 0), mgp = c(2, 1, 0))
plot(coh[,1]~freqs, type='l', 
     ylab = 'WL and Barometric',
     xlab = 'Frequency (cpd)',
     xlim = rng_ba,
     main = 'Coherence')
points(y = coh_base[,1], x = freqs, type='l',  col = 'red')
legend('bottomleft', 
       legend = c('psd', 'base'), 
       col = c('black', 'red'), lty = 1)

## ----eval=TRUE, echo=FALSE, label="coh_et_plot", fig.height=3.3---------------
par(mar = c(3, 3, 2, 0), mgp = c(2, 1, 0))

plot(coh[,2]~freqs, type='l', 
     ylab = 'WL and Earth tide',
     xlab = 'Frequency (cpd)',
     xlim = rng_et,
     main = 'Coherence')
points(y = coh_base[,2], x = freqs, type='l',  col = 'red')



## ----eval=TRUE, echo=TRUE, label="phase"--------------------------------------
phase <- pspectrum(dat, riedsid_column= 0L, verbose = FALSE,
                 niter = 0, ntap.init = 11)$phase
phase_base <- spec.pgram(dat[, c(3,1,2)], spans = c(11), fast = FALSE, plot = FALSE)$phase

## ----eval=TRUE, echo=FALSE, label="phase_ba_plot", fig.height=3.3-------------
par(mar = c(3, 3, 2, 0), mgp = c(2, 1, 0))
plot(phase[,1]~freqs, type='l', 
     ylab = 'WL and Barometric',
     xlab = 'Frequency (cpd)',
     xlim = rng_ba,
     main = 'Phase Difference')
points(y = phase_base[,1], x = freqs, type='l',  col = 'red')
legend('bottomleft', 
       legend = c('psd', 'base'), 
       col = c('black', 'red'), lty = 1)

## ----eval=TRUE, echo=FALSE, label="phase_et_plot", fig.height=3.3-------------
par(mar = c(3, 3, 2, 0), mgp = c(2, 1, 0))

plot(phase[,2]~freqs, type='l', 
     ylab = 'WL and Earth tide',
     xlab = 'Frequency (cpd)',
     xlim = rng_et,
     main = 'Phase')
points(y = phase_base[,2], x = freqs, type='l',  col = 'red')



## ----eval=TRUE, echo=TRUE, label="transfer gain", fig.height=3.3--------------
transfer <- pspectrum(dat, riedsid_column= -1L, verbose = FALSE)$transfer
gain <- Mod(transfer)

## ----eval=TRUE, echo=FALSE, label="transfer gain plot", fig.height=3.3--------


par(mar = c(3, 3, 2, 0), mgp = c(2, 1, 0))
plot(gain[, 1]~freqs, type='l', 
     xlim = rng_ba, 
     ylim = c(0, 1),
     ylab = 'WL and Barometric',
     xlab = 'Frequency (cpd)',
     main = 'Gain')
plot(gain[, 2]~freqs, type='l', 
     xlim = rng_et,
     ylim = c(0, 0.00005),
     ylab = 'WL and Earth tide',
     xlab = 'Frequency (cpd)',
     main = 'Gain')

## ----eval=TRUE, echo=TRUE, label="transfer phase", fig.height=3.3-------------
phase <- Arg(transfer)

## ----eval=TRUE, echo=FALSE, label="transfer phase plot", fig.height=3.3-------
par(mar = c(3, 3, 2, 0), mgp = c(2, 1, 0))

plot(phase[, 1]~freqs, type='l', xlim = rng_ba,
     ylab = 'WL and Barometric',
     xlab = 'Frequency (cpd)',
     main = 'Phase')
plot(phase[, 2]~freqs, type='l', xlim = rng_et, 
     ylab = 'WL and Earth tide',
     xlab = 'Frequency (cpd)',
     main = 'Phase')

## ----eval=TRUE, echo=TRUE, label=SI-------------------------------------------
utils::sessionInfo()

