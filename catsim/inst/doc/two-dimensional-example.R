## ----init, include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = "styler"
)

## ----setup--------------------------------------------------------------------
library(catsim)

## ----besag, fig.width = 6, fig.height = 2-------------------------------------
set.seed(20200308)

multfactor <- 3
onesmat <- matrix(1, multfactor, multfactor)
exmat <- (besag %% 2) ### converting to binary - it's ones and twos

bigmat <- exmat %x% onesmat

bigmat <- exmat %x% onesmat
shift <- 6
basemat <- bigmat[1:(96 * multfactor), 1:(84 * multfactor)]
shiftmat <- bigmat[(1 + shift):(96 * multfactor + shift), 1:(84 * multfactor)]
acc <- mean(basemat == shiftmat)
errmat <- (matrix(
  sample(0:1, 96 * 84 * (multfactor^2),
    replace = TRUE, prob = c(acc, ((1 - acc)))
  ),
  (96 * multfactor), (84 * multfactor)
) + basemat) %% 2
par(mfrow = c(1,3), mar = c(0,0,0,0))
image(basemat, xaxt = "n", yaxt = "n"); grid(lwd = 3, nx = 6, col = "black")
image(shiftmat, xaxt = "n", yaxt = "n"); grid(lwd = 3, nx = 6, col = "black")
image(errmat, xaxt = "n", yaxt = "n"); grid(lwd = 3, nx = 6, col = "black")



## ----besagrating--------------------------------------------------------------
catsim(basemat, shiftmat)

catsim(basemat, errmat)


## ----besagjacc----------------------------------------------------------------
catsim(basemat, shiftmat, method = "Jaccard")

catsim(basemat, errmat, method = "Jaccard")


## ----construction, fig.width = 2, fig.height = 2------------------------------

multfactor <- 20
onesmat <- matrix(1, multfactor, multfactor)
exmat <- matrix(0, 12, 12)
exmat[2:10, c(2, 3, 9, 10)] <- 1
exmat[2:10, c(4, 5, 7, 8)] <- 2
exmat[5:7, 2:10] <- 1
exmat[2:10, c(6)] <- 3
par( mar = c(0,0,0,0))
image(exmat[11:1, 1:11], xaxt = "n", yaxt = "n")
  grid(col = "black", lwd = 3)
bigmat <- exmat %x% onesmat


## ----distortion, fig.height = 2, fig.width = 4--------------------------------
set.seed(20200323)
shift  <-  6 # we are shifting horizontally by six pixels
basemat <- bigmat[1:(11 * multfactor), 1:(11 * multfactor)]
shiftmat <- bigmat[(1 + shift) : (11 * multfactor + shift), 1: (11 * multfactor)]
### Here we have shifted the matrix slightly
par(mfrow = c(1, 2), mar = c(0,0,0,0))
image(shiftmat, xaxt = "n", yaxt = "n")
grid(col = "black", lwd = 3)

### computing the error rate
acc  <-  mean(basemat == shiftmat)
errmat <- (matrix(sample(0:3, 121 * (multfactor^2), replace = TRUE,
                         prob = c(acc, rep((1 - acc) / 3, 3))),
                        (11 * multfactor), (11 * multfactor)) + basemat) %% 4

### here we have made an image that matches its accuracy with salt and pepper noise
image(errmat, xaxt = "n", yaxt = "n")
grid(col = "black", lwd = 3)


## ----onelevel-----------------------------------------------------------------
library(catsim)
### comparing base to shifted
catsim(basemat, shiftmat, weights = 1)

### comparing base to salt-and-pepper
catsim(basemat, errmat, weights = 1)

### looking at accuracy
mean(basemat == shiftmat)
mean(basemat == errmat)


## ----fivelevels---------------------------------------------------------------
catsim(basemat, shiftmat, weights = rep(.2, 5))

### comparing base to salt-and-pepper
catsim(basemat, errmat, weights = rep(.2, 5))


## ----differentmetrics---------------------------------------------------------
catsim(basemat, shiftmat, weights = rep(.5,2), method="Rand")
catsim(basemat, shiftmat, weights = rep(.5,2), method="NMI")


## ----windowsize---------------------------------------------------------------
catsim(basemat, shiftmat, weights = rep(1 / 3, 3), window = 20)
catsim(basemat, shiftmat, weights = rep(1 / 3, 3), window = 5)
catsim(basemat, errmat, weights = rep(1 / 3, 3), window = 20)
catsim(basemat, errmat, weights = rep(1 / 3, 3), window = 5)


## ----windowarg----------------------------------------------------------------
catsim(basemat, shiftmat, weights = rep(.2, 5), window = c(11, 5))

