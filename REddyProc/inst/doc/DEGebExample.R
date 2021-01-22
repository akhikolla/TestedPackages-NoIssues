## ---- include = FALSE----------------------------------------------------
# do not execute on CRAN: 
# https://stackoverflow.com/questions/28961431/computationally-heavy-r-vignettes
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_",
             "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
knitr::opts_chunk$set(eval = !is_check)

## ----setup, include = FALSE----------------------------------------------
library(knitr)
# rmarkdown::render("vignettes/DEGebExample.Rmd")
opts_chunk$set(out.extra = 'style = "display:block; margin: auto"'
    #, fig.align = "center"
    , fig.width = 4.3, fig.height = 3.2, dev.args = list(pointsize = 10)
    , message = FALSE
    )
knit_hooks$set(spar = function(before, options, envir) {
    if (before) {
        par( las = 1 )                   #also y axis labels horizontal
        par(mar = c(2.0,3.3,0,0) + 0.3 )  #margins
        par(tck = 0.02 )                          #axe-tick length inside plots             
        par(mgp = c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
     }
})
# genVigs("DEGebExample")

## ----results = 'hide'----------------------------------------------------
#isDevelopMode <- TRUE
if (!exists("isDevelopMode")) library(REddyProc)
set.seed(0815)      # for reproducible results

## ------------------------------------------------------------------------
data(DEGebExample)
summary(DEGebExample)

## ------------------------------------------------------------------------
DEGebExample$VPD <- fCalcVPDfromRHandTair(DEGebExample$rH, DEGebExample$Tair)
EProc <- sEddyProc$new('DE-Geb', DEGebExample, c('NEE','Rg','Tair','VPD', 'Ustar'))
EProc$sSetLocationInfo(LatDeg = 51.1, LongDeg = 10.9, TimeZoneHour = 1)  #Location of Gebesee

## ----DEGeb_estUStar1a, spar = TRUE, fig.width = 10-----------------------
seasonStarts <- as.data.frame( do.call( rbind, list(
	  c(70,2004)
          ,c(210,2004)
          ,c(320,2004)
          ,c(70,2005)
          ,c(180,2005)
          ,c(320,2005)
          ,c(120,2006)
          ,c(305,2006) 		
)))
seasonFactor <- usCreateSeasonFactorYdayYear(
  DEGebExample$DateTime - 15*60, starts = seasonStarts)
plot( NEE ~ DateTime, DEGebExample )
seasonStartsDate <- fConvertTimeToPosix( data.frame(Year = seasonStarts[,2]
	, DoY = seasonStarts[,1], Hour = 0.25), 'YDH'
	, Year = "Year", Day = "DoY", Hour = "Hour")
abline( v = seasonStartsDate$DateTime)

## ----DEGeb_estUStar1b, cache = TRUE--------------------------------------
(uStarTh <- EProc$sEstUstarThold(seasonFactor = seasonFactor))
#EProc$useSeaonsalUStarThresholds()
# estimation can be inspected by plotting the saturation of NEE with UStar 
# for the temperatures of one season
#EProc$sPlotNEEVersusUStarForSeason( levels(uStarTh$season)[2] )

## ----DEGeb_estUStar1c----------------------------------------------------
EProc$useSeaonsalUStarThresholds()
EProc$sGetUstarScenarios()

## ----DEGeb_gapFill1, cache = TRUE, eval = FALSE, output = 'hide', message = FALSE----
#  EProc$sMDSGapFillAfterUstar('NEE', FillAll = FALSE, isVerbose = FALSE)

## ----DEGeb_estUStarBoot1, cache = TRUE, message = FALSE------------------
EProc <- sEddyProc$new('DE-Geb', DEGebExample, c('NEE','Rg','Tair','VPD', 'Ustar'))
EProc$sSetLocationInfo(LatDeg = 51.1, LongDeg = 10.9, TimeZoneHour = 1)  #Location of Gebesee
# here, because of processing time only 30 samples instead of 100, and 10% and 90% 
# percentile instead of default 5%,50%, and 95% with 100 samples
EProc$sEstimateUstarScenarios( 
  seasonFactor = seasonFactor, nSample = 30L, probs = c(0.1,0.9))
#(uStarScens <- usGetSeasonalSeasonUStarMap(
#  EProc$sGetEstimatedUstarThresholdDistribution()
#))
#EProc$sSetUstarScenarios(uStarScens)
EProc$useSeaonsalUStarThresholds()
EProc$sGetUstarScenarios()

## ----DEGeb_gapFillBoot1, cache = TRUE, message = FALSE, output = 'hide'----
EProc$sMDSGapFillUStarScens('NEE', FillAll = FALSE)

## ------------------------------------------------------------------------
grep("^NEE.*_f$", colnames( EProc$sExportResults()), value = TRUE )

## ----DEGeb_fluxPart1, cache = TRUE, message = FALSE----------------------
EProc$sMDSGapFill('Tair', FillAll = FALSE)
EProc$sApplyUStarScen( EProc$sMRFluxPartition )
#grep("U10", colnames(EProc$sExportResults()), value = TRUE) 	
grep("^GPP.*_f$", colnames( EProc$sExportResults()), value = TRUE )

## ----DEGeb_estUStarCPT, cache = TRUE-------------------------------------
EProc <- sEddyProc$new(
  'DE-Geb', DEGebExample, c('NEE','Rg','Tair','VPD', 'Ustar'))
resUStar <- EProc$sEstUstarThold(
					ctrlUstarEst = usControlUstarEst(isUsingCPTSeveralT = TRUE)
					, seasonFactor = seasonFactor
			)
#(uStarThCP <- usGetSeasonalSeasonUStarMap(resUStar))
EProc$useSeaonsalUStarThresholds()
EProc$sGetUstarScenarios()

