## ---- include = FALSE----------------------------------------------------
# do not execute on CRAN: 
# https://stackoverflow.com/questions/28961431/computationally-heavy-r-vignettes
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_",
             "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
knitr::opts_chunk$set(eval = !is_check)

## ----setup, include = FALSE----------------------------------------------
library(knitr)
#rmarkdown::render("vignettes/uStarCases.Rmd","md_document")
opts_knit$set(root.dir = '..')
opts_chunk$set(
    #, fig.align = "center"
    #, fig.width = 3.27, fig.height = 2.5, dev.args = list(pointsize = 10)
    #,cache = TRUE
    #, fig.width = 4.3, fig.height = 3.2, dev.args = list(pointsize = 10)
    #, fig.width = 6.3, fig.height = 6.2, dev.args = list(pointsize = 10)
    # works with html but causes problems with latex
    #,out.extra = 'style = "display:block; margin: auto"' 
    )
knit_hooks$set(spar = function(before, options, envir) {
    if (before) {
        par(las = 1 )                   #also y axis labels horizontal
        par(mar = c(2.0,3.3,0,0) + 0.3 )  #margins
        par(tck = 0.02 )                          #axe-tick length inside plots             
        par(mgp = c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
     }
})

## ---- include = FALSE----------------------------------------------------
#themeTw <- theme_bw(base_size = 10) + theme(axis.title = element_text(size = 9))

## ----init, message = FALSE, output = 'hide'------------------------------
#+++ load libraries used in this vignette
library(REddyProc)
library(dplyr)
#+++ define directory for outputs
outDir <- tempdir()  # CRAN policy dictates to write only to this dir in examples
#outDir <- "out"     # to write to subdirectory of current users dir
#+++ Add time stamp in POSIX time format to example data 
# and filter long runs of equal NEE values
EddyDataWithPosix <- fConvertTimeToPosix(
  filterLongRuns(Example_DETha98, "NEE")
  , 'YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')

## ----noUStar, message = FALSE--------------------------------------------
EProc <- sEddyProc$new(
  'DE-Tha', EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
EProc$sMDSGapFill('NEE')
grep("NEE.*_f$",names(EProc$sExportResults()), value = TRUE)

## ----userUStar, message = FALSE------------------------------------------
EProc <- sEddyProc$new(
  'DE-Tha', EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
uStar <- 0.46
EProc$sMDSGapFillAfterUstar('NEE', uStarTh = uStar)
grep("NEE.*_f$",names(EProc$sExportResults()), value = TRUE)

## ----singleUStar, message = FALSE----------------------------------------
EProc <- sEddyProc$new(
  'DE-Tha', EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
# estimating the thresholds based on the data (without bootstrap)
(uStarTh <- EProc$sEstUstarThold())
# may plot saturation of NEE with UStar for a specified season to pdf
EProc$sPlotNEEVersusUStarForSeason(levels(uStarTh$season)[3], dir = outDir )

## ----singleUStarGapfill, message = FALSE---------------------------------
#EProc$useAnnualUStarThresholds()
EProc$sMDSGapFillAfterUstar('NEE')
grep("NEE.*_f$",names(EProc$sExportResults()), value = TRUE)

## ----uStarScen, results='hold'-------------------------------------------
EProc <- sEddyProc$new(
  'DE-Tha', EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
EProc$sEstimateUstarScenarios(
    nSample = 100L, probs = c(0.05, 0.5, 0.95))
# inspect the thresholds to be used by default
EProc$sGetUstarScenarios()

## ------------------------------------------------------------------------
(uStarThAgg <- EProc$sGetEstimatedUstarThresholdDistribution())

## ----uStarScenSetSeasonal------------------------------------------------
#EProc$sSetUstarScenarios(
#  usGetSeasonalSeasonUStarMap(uStarThAgg)[,-2])
EProc$useSeaonsalUStarThresholds()
# inspect the changed thresholds to be used
EProc$sGetUstarScenarios()

## ----uStarScenGapfill, message=FALSE-------------------------------------
EProc$sMDSGapFillUStarScens("NEE")
grep("NEE_.*_f$",names(EProc$sExportResults()), value = TRUE)

## ----uStarScenMRPart, message=FALSE--------------------------------------
EProc$sSetLocationInfo(LatDeg = 51.0, LongDeg = 13.6, TimeZoneHour = 1)
EProc$sMDSGapFill('Tair', FillAll = FALSE, minNWarnRunLength = NA)
EProc$sMDSGapFill('Rg', FillAll = FALSE, minNWarnRunLength = NA)
EProc$sMDSGapFill('VPD', FillAll = FALSE, minNWarnRunLength = NA)
EProc$sMRFluxPartitionUStarScens()
grep("GPP_.*_f$",names(EProc$sExportResults()), value = TRUE)
if (FALSE) {
  # run only interactively, because it takes long
  EProc$sGLFluxPartitionUStarScens(uStarScenKeep = "U50")
  grep("GPP_DT_.*_f$",names(EProc$sExportResults()), value = TRUE)
}

