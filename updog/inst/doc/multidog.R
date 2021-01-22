## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=4.5,
  fig.height=3.5
)

## ----setup--------------------------------------------------------------------
library(updog)
data("uitdewilligen")

## -----------------------------------------------------------------------------
refmat  <- t(uitdewilligen$refmat)
sizemat <- t(uitdewilligen$sizemat)
ploidy  <- uitdewilligen$ploidy

## -----------------------------------------------------------------------------
setdiff(colnames(sizemat), colnames(refmat))
setdiff(rownames(sizemat), rownames(refmat))

## -----------------------------------------------------------------------------
parallel::detectCores()

## -----------------------------------------------------------------------------
mout <- multidog(refmat = refmat, 
                 sizemat = sizemat, 
                 ploidy = ploidy, 
                 model = "norm",
                 nc = 2)

## -----------------------------------------------------------------------------
plot(mout, indices = c(1, 5, 100))

## -----------------------------------------------------------------------------
str(mout$snpdf)

## -----------------------------------------------------------------------------
str(mout$inddf)

## -----------------------------------------------------------------------------
genomat <- format_multidog(mout, varname = "geno")
head(genomat)

## -----------------------------------------------------------------------------
dim(mout$snpdf)
dim(mout$inddf)
mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > exp(-1) & bias < exp(1))
dim(mout_cleaned$snpdf)
dim(mout_cleaned$inddf)

