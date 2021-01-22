## ----include=FALSE------------------------------------------------------------
library(knitr)
opts_chunk$set(
tidy=FALSE,
size="small"
)

## ----echo=FALSE---------------------------------------------------------------
#require( gRbase )
prettyVersion <- packageDescription("gRbase")$Version
prettyDate <- format(Sys.Date())

## ----echo=F---------------------------------------------------------------------------------------
library(gRbase)
options("width"=100, "digits"=4)
options(useFancyQuotes="UTF-8")
#chk = 'markup'
chk = "hide"

## -------------------------------------------------------------------------------------------------
hec <- c(32, 53, 11, 50, 10, 25, 36, 66, 9, 34, 5, 29) 
dim(hec) <- c(2, 3, 2)
dimnames(hec) <- list(Hair = c("Black", "Brown"), 
                      Eye = c("Brown", "Blue", "Hazel"), 
                      Sex = c("Male", "Female"))
hec

## -------------------------------------------------------------------------------------------------
##flat <- function(x) {ftable(x, row.vars=1)}
flat <- function(x, n=4) {as.data.frame.table(x) %>% head(n)}
hec %>% flat

## -------------------------------------------------------------------------------------------------
is.named.array(hec)

## -------------------------------------------------------------------------------------------------
dn <- list(Hair=c("Black", "Brown"), Eye=~Brown:Blue:Hazel, Sex=~Male:Female)
counts <- c(32, 53, 11, 50, 10, 25, 36, 66, 9, 34, 5, 29)
z3 <- tabNew(~Hair:Eye:Sex, levels=dn, value=counts) 
z4 <- tabNew(c("Hair", "Eye", "Sex"), levels=dn, values=counts)

## -------------------------------------------------------------------------------------------------
z5 <- tabNew(~Hair:Eye:Sex, levels=c(2, 3, 2), values = counts)
dimnames(z5) %>% str

## -------------------------------------------------------------------------------------------------
z6 <- tabNew(~Hair:Eye:Sex, levels=c(2, 3, 2), values=counts, normalize="first")
z6 %>% flat

## -------------------------------------------------------------------------------------------------
tabNormalize(z5, "first") %>% flat

## ----results=chk----------------------------------------------------------------------------------
tabSlice(hec, slice=list(Eye=c("Blue", "Hazel"), Sex="Female"))
## Notice: levels can be written as numerics
## tabSlice(hec, slice=list(Eye=2:3, Sex="Female"))

## -------------------------------------------------------------------------------------------------
tabSlice(hec, slice=list(Eye=c("Blue", "Hazel"), Sex="Female"), drop=FALSE)

## ----results=chk----------------------------------------------------------------------------------
## A vector:
t1 <- tabSlice(hec, slice=list(Hair=1, Sex="Female")); t1
## A 1-dimensional array:
t2 <- tabSlice(hec, slice=list(Hair=1, Sex="Female"), as.array=TRUE); t2 
## A higher dimensional array (in which some dimensions only have one level)
t3 <- tabSlice(hec, slice=list(Hair=1, Sex="Female"), drop=FALSE); t3

## -------------------------------------------------------------------------------------------------
t2 %>% flat
t3 %>% flat

## -------------------------------------------------------------------------------------------------
he <- tabMarg(hec, c("Hair", "Eye"))
he

## ----results=chk----------------------------------------------------------------------------------
## Alternatives
tabMarg(hec, ~Hair:Eye)
tabMarg(hec, c(1, 2))
hec %a_% ~Hair:Eye

## -------------------------------------------------------------------------------------------------
he1 <- tabMarg(hec, c("Hair", "Eye"))
he2 <- tabMarg(he1, c("Hair", "Eye"))
tabEqual(he1, he2)

## -------------------------------------------------------------------------------------------------
extra.dim <- list(Sex=c("Male", "Female"))
tabExpand(he, extra.dim) 

## ----results=chk----------------------------------------------------------------------------------
## Alternatives
he %a^% extra.dim

## -------------------------------------------------------------------------------------------------
(he %a^% extra.dim) %a_% c("Hair", "Eye")

## -------------------------------------------------------------------------------------------------
tabPerm(hec, ~Eye:Sex:Hair) %>% flat 

## ----results=chk----------------------------------------------------------------------------------
tabPerm(hec, c("Eye", "Sex", "Hair"))
tabPerm(hec, c(2, 3, 1)) 
tabPerm(hec, ~Ey:Se:Ha) 
tabPerm(hec, c("Ey", "Se", "Ha"))

## -------------------------------------------------------------------------------------------------
hec2 <- tabPerm(hec, 3:1)
tabEqual(hec, hec2)

## ----results=chk----------------------------------------------------------------------------------
## Alternative
hec %a==% hec2

## ----results=chk----------------------------------------------------------------------------------
hec2 <- tabPerm(hec, 3:1)
tabAlign(hec2, hec)

## ----result=chk-----------------------------------------------------------------------------------
## Alternative:
tabAlign(hec2, dimnames(hec))

## -------------------------------------------------------------------------------------------------
hs <- tabMarg(hec, ~Hair:Eye)
tabMult(he, hs)

## ----results=chk----------------------------------------------------------------------------------
tabAdd(he, hs) 
tabSubt(he, hs)
tabMult(he, hs)
tabDiv(he, hs) 
tabDiv0(he, hs) ## Convention 0/0 = 0

## ----results=chk----------------------------------------------------------------------------------
## Alternative
he %a+% hs
he %a-% hs
he %a*% hs
he %a/% hs
he %a/0% hs ## Convention 0/0 = 0

## -------------------------------------------------------------------------------------------------
es <- tabMarg(hec, ~Eye:Sex)
tabSum(he, hs, es)  
## tabSum(list(he, hs, es))

## -------------------------------------------------------------------------------------------------
tabDist(hec, marg=~Hair:Eye)
tabDist(hec, cond=~Sex) 
tabDist(hec, marg=~Hair, cond=~Sex) 

## -------------------------------------------------------------------------------------------------
tabSliceMult(es, list(Sex="Female"), val=10, comp=0)

## -------------------------------------------------------------------------------------------------
yn <- c("y","n")
lev <- list(rain=yn, sprinkler=yn, wet=yn)
r <- tabNew(~rain, levels=lev, values=c(.2, .8))
s_r <- tabNew(~sprinkler:rain, levels = lev, values = c(.01, .99, .4, .6))
w_sr <- tabNew( ~wet:sprinkler:rain, levels=lev, 
             values=c(.99, .01, .8, .2, .9, .1, 0, 1))
r 
s_r  %>% flat
w_sr %>% flat

## -------------------------------------------------------------------------------------------------
joint <- tabProd(r, s_r, w_sr); joint %>% flat

## -------------------------------------------------------------------------------------------------
tabDist(joint, marg=~rain, cond=~wet)

## ----results='hide'-------------------------------------------------------------------------------
## Alternative:
rw <- tabMarg(joint, ~rain + wet)
tabDiv(rw, tabMarg(rw, ~wet))
## or
rw %a/% (rw %a_% ~wet)

## -------------------------------------------------------------------------------------------------
## Alternative:
x <- tabSliceMult(rw, slice=list(wet="y")); x
tabDist(x, marg=~rain)

## -------------------------------------------------------------------------------------------------
data(lizard, package="gRbase")
lizard %>% flat

## -------------------------------------------------------------------------------------------------
myips <- function(indata, glist){
    fit   <- indata
    fit[] <-  1
    ## List of sufficient marginal tables
    md    <- lapply(glist, function(g) tabMarg(indata, g))

    for (i in 1:4){
        for (j in seq_along(glist)){
            mf  <- tabMarg(fit, glist[[j]])
            # adj <- tabDiv( md[[ j ]], mf)
            # fit <- tabMult( fit, adj )
            ## or
            adj <- md[[ j ]] %a/% mf
            fit <- fit %a*% adj
        }
    }
    pearson <- sum((fit - indata)^2 / fit)
    list(pearson=pearson, fit=fit)
}

glist <- list(c("species", "diam"),c("species", "height"),c("diam", "height"))

fm1 <- myips(lizard, glist)
fm1$pearson
fm1$fit %>% flat

fm2 <- loglin(lizard, glist, fit=T)
fm2$pearson
fm2$fit %>% flat

## -------------------------------------------------------------------------------------------------
hec
c(hec)

## -------------------------------------------------------------------------------------------------
cell2name <- function(cell, dimnames){
    unlist(lapply(1:length(cell), function(m) dimnames[[m]][cell[m]]))
}
cell2name(c(2,3,1), dimnames(hec))

## -------------------------------------------------------------------------------------------------
cell2entry(c(2,3,1), dim=c(2, 3, 2))
entry2cell(6, dim=c(2, 3, 2))

## -------------------------------------------------------------------------------------------------
next_cell(c(2,3,1), dim=c(2, 3, 2))

## -------------------------------------------------------------------------------------------------
next_cell_slice(c(1,3,1), slice_marg=2, dim=c( 2, 3, 2 ))
next_cell_slice(c(2,3,1), slice_marg=2, dim=c( 2, 3, 2 ))

## -------------------------------------------------------------------------------------------------
slice2entry(slice_cell=3, slice_marg=2, dim=c( 2, 3, 2 ))

## -------------------------------------------------------------------------------------------------
r <- slice2entry(slice_cell=3, slice_marg=2, dim=c( 2, 3, 2 ))
lapply(lapply(r, entry2cell, c( 2, 3, 2 )),
       cell2name, dimnames(hec))

## -------------------------------------------------------------------------------------------------
head( fact_grid( c(2, 3, 2) ), 6 )

## -------------------------------------------------------------------------------------------------
head( expand.grid(list(1:2, 1:3, 1:2)), 6 )

## -------------------------------------------------------------------------------------------------
hec[, 2:3, ]  %>% flat  ## A 2 x 2 x 2 array
hec[1, , 1]             ## A vector
hec[1, , 1, drop=FALSE] ## A 1 x 3 x 1 array

## ----results=chk----------------------------------------------------------------------------------
do.call("[", c(list(hec), list(TRUE, 2:3, TRUE)))  %>% flat
do.call("[", c(list(hec), list(1, TRUE, 1))) 
do.call("[", c(list(hec), list(1, TRUE, 1), drop=FALSE)) 

## ----results=chk----------------------------------------------------------------------------------
tabSlicePrim(hec, slice=list(TRUE, 2:3, TRUE))  %>% flat
tabSlice(hec, slice=list(c(2, 3)), margin=2) %>% flat

tabSlicePrim(hec, slice=list(1, TRUE, 1))  
tabSlice(hec, slice=list(1, 1), margin=c(1, 3)) 

tabSlicePrim(hec, slice=list(1, TRUE, 1), drop=FALSE)  
tabSlice(hec, slice=list(1, 1), margin=c(1, 3), drop=FALSE) 

