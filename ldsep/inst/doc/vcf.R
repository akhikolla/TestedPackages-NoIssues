## ---- include = FALSE---------------------------------------------------------
runthis <- requireNamespace("VariantAnnotation", quietly = TRUE) &&
  requireNamespace("updog", quietly = TRUE)

if (runthis) {
  runthis <- packageVersion("updog") >= "2.0.2"
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 4,
  eval = runthis
)
set.seed(1)

## ----setup, message=FALSE-----------------------------------------------------
library(VariantAnnotation)
library(updog)
library(ldsep)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages(c("ldsep", "updog", "BiocManager"))
#  BiocManager::install("VariantAnnotation")

## -----------------------------------------------------------------------------
packageVersion("updog")

## -----------------------------------------------------------------------------
uit_file <- system.file("extdata", "subuit.vcf", 
                        package = "ldsep",
                        mustWork = TRUE)

## -----------------------------------------------------------------------------
subuit <- readVcf(file = uit_file)
class(subuit)

## -----------------------------------------------------------------------------
geno(header(subuit))

## -----------------------------------------------------------------------------
sizemat <- geno(subuit)$DP
refmat <- geno(subuit)$RA

## -----------------------------------------------------------------------------
ploidy <- 4
mout <- multidog(refmat = refmat, 
                 sizemat = sizemat, 
                 ploidy = ploidy, 
                 model = "norm")

## ---- warning=FALSE, fig.show="hold", out.width="25%", results='hide'---------
plot(mout, indices = sample(1:nrow(subuit), 3))

## -----------------------------------------------------------------------------
msub <- filter_snp(x = mout, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95)
nrow(msub$snpdf)

## -----------------------------------------------------------------------------
varnames <- paste0("logL_", 0:ploidy)
varnames
larray <- format_multidog(x = msub, varname = varnames)
class(larray)
dim(larray)

## -----------------------------------------------------------------------------
pmmat <- format_multidog(x = msub, varname = "postmean")
class(pmmat)
dim(pmmat)

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages(c("fitPoly", "reshape2"))
#  library(fitPoly)
#  library(reshape2)

## ---- eval = FALSE------------------------------------------------------------
#  ratiodf <- melt(data = refmat / sizemat,
#                  varnames = c("MarkerName", "SampleName"),
#                  value.name = "ratio")
#  saveMarkerModels(ploidy = ploidy,
#                   data = ratiodf,
#                   filePrefix = "uit",
#                   rdaFiles = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  load("./uit_scores.RData")
#  load("./uit_models.RData")

## ----eval=FALSE---------------------------------------------------------------
#  modeldata <- subset(x = modeldata, subset = pmax(P0, P1, P2, P3, P4) < 0.95)
#  scores <- scores[scores$MarkerName %in% modeldata$MarkerName, ]

## ---- eval=FALSE--------------------------------------------------------------
#  mergedf <- merge(x = scores[c("MarkerName", "SampleName", paste0("P", 0:ploidy))],
#                   y = modeldata[c("MarkerName", paste0("P", 0:ploidy))],
#                   by = "MarkerName",
#                   all.x = TRUE)
#  for (i in 0:ploidy) {
#    current_post <- paste0("P", i, ".x")
#    current_prior <- paste0("P", i, ".y")
#    current_ll <- paste0("logL_", i)
#    mergedf[[current_ll]] <- log(mergedf[[current_post]]) - log(mergedf[[current_prior]])
#  }
#  
#  fp_llarray <- acast(
#    melt(data = mergedf[c("MarkerName", "SampleName", paste0("logL_", 0:ploidy))],
#         id.vars = c("MarkerName", "SampleName"),
#         variable.name = "dosage",
#         value.name = "logl"),
#    formula = "MarkerName ~ SampleName ~ dosage",
#    value.var = "logl",
#  )

## -----------------------------------------------------------------------------
like_ld <- mldest(geno = larray, K = ploidy, type = "comp")

## ---- out.width="50%"---------------------------------------------------------
plot(like_ld)

## -----------------------------------------------------------------------------
mom_ld <- mldest(geno = pmmat, K = ploidy, type = "comp")

## ---- out.width="50%"---------------------------------------------------------
plot(mom_ld)

## -----------------------------------------------------------------------------
par(mar = c(2.4, 2.8, 0, 0) + 0.5, mgp = c(1.8, 0.6, 0))
plot(mom_ld$r2, like_ld$r2, 
     xlab = expression(paste(textstyle(Naive), ~~hat(r)^2)), 
     ylab = expression(paste(textstyle(MLE), ~~hat(r)^2)), 
     pch  = 20)
abline(0, 1, lty = 2, col = 2)

## -----------------------------------------------------------------------------
ldmat <- format_lddf(obj = like_ld, element = "r2")
ldmat[1:4, 1:4]

