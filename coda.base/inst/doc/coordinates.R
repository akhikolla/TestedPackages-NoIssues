## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(coda.base)
# By default basis is not shown, in this vignette we turn on basis showing.
options('coda.base.basis' = TRUE)
data('parliament2017')
X = parliament2017[,c('erc','jxcat','psc','cs')]

## -----------------------------------------------------------------------------
H1.alr = coordinates(X, basis = 'alr')
head(H1.alr)

## -----------------------------------------------------------------------------
alr_basis(dim = 4)

## -----------------------------------------------------------------------------
B.alr = alr_basis(dim = 4, numerator = c(4,2,3), denominator = 1)
B.alr

## -----------------------------------------------------------------------------
H2.alr = coordinates(X, basis = B.alr)
head(H2.alr)

## -----------------------------------------------------------------------------
H.clr = coordinates(X, basis = 'clr')
head(H.clr)

## -----------------------------------------------------------------------------
H1.ilr = coordinates(X)
head(H1.ilr)

## -----------------------------------------------------------------------------
all.equal( coordinates(X, basis = 'ilr'),
           H1.ilr )

## ---- fig.width=5.5, fig.height=4, fig.align='center', caption='Variance of principal components coordinates'----
H2.ilr = coordinates(X, basis = 'pc')
head(H2.ilr)
barplot(apply(H2.ilr, 2, var))

## -----------------------------------------------------------------------------
cov(H2.ilr)

## ---- fig.width=5.5, fig.height=4, fig.align='center', caption='Variance of principal balances coordinates'----
H3.ilr = coordinates(X, basis = 'pb')
head(H3.ilr)
barplot(apply(H3.ilr, 2, var))

## -----------------------------------------------------------------------------
cor(H3.ilr)

## -----------------------------------------------------------------------------
X100 = exp(matrix(rnorm(1000*100), ncol = 100))

## -----------------------------------------------------------------------------
PB1.ward = pb_basis(X100, method = 'cluster')

## -----------------------------------------------------------------------------
PB1.constrained = pb_basis(X100, method = 'constrained')

## -----------------------------------------------------------------------------
PC_approx = coordinates(X100, cbind(pc_basis(X100)[,1], PB1.ward[,1], PB1.constrained[,1]))
names(PC_approx) = c('PC', 'Ward', 'Constrained')
apply(PC_approx, 2, var)

## -----------------------------------------------------------------------------
H4.ilr = coordinates(X, basis = 'cdp')
head(H4.ilr)

## -----------------------------------------------------------------------------
B = matrix(c(-1,-1,2,0,
             1,0,-0.5,-0.5,
             -0.5,0.5,0,0), ncol = 3)
H1.man = coordinates(X, basis = B)
head(H1.man)

## ---- eval=FALSE--------------------------------------------------------------
#  B.man = sbp_basis(b1 = erc~jxcat,
#                    b2 = psc~cs,
#                    b3 = erc+jxcat~psc+cs,
#                    data=X)
#  H2.man = coordinates(X, basis = B.man)
#  head(H2.man)

## -----------------------------------------------------------------------------
B = sbp_basis(b1 = erc+jxcat~psc+cs, 
              data=X)
H3.man = coordinates(X, basis = B)
head(H3.man)

## -----------------------------------------------------------------------------
B = sbp_basis(b1 = erc~jxcat+psc~cs, 
              b2 = jxcat~erc+psc+cs,
              b3 = psc~erc+jxcat+cs,
              b4 = cs~erc+jxcat+psc,
              data=X)
H4.man = coordinates(X, basis = B)
head(H4.man)

## -----------------------------------------------------------------------------
P =  matrix(c(1, 1,-1,-1,
              1,-1, 0, 0,
              0, 0, 1,-1), ncol= 3)
B = sbp_basis(P)
H5.man = coordinates(X, basis = B)
head(H5.man)

