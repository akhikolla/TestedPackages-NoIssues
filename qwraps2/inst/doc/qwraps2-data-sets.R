## -----------------------------------------------------------------------------
library(qwraps2)
packageVersion("qwraps2")

## -----------------------------------------------------------------------------
mtcars2 <- mtcars
str(mtcars2)

## -----------------------------------------------------------------------------
mtcars2$cyl_character <- paste(mtcars2$cyl, "cylinders")
mtcars2$cyl_factor    <- factor(mtcars2$cyl,
                                levels = c(6, 4, 8),
                                labels = paste( c(6, 4, 8), "cylinders"))

## -----------------------------------------------------------------------------
mtcars2$gear_factor <-
  factor(mtcars2$gear, levels = c(3, 4, 5), labels = paste(c(3, 4, 5), "forward gears"))

## -----------------------------------------------------------------------------
mtcars2$engine <-
  factor(mtcars2$vs, levels = c(0, 1), labels = c("V-shaped", "straight"))

## -----------------------------------------------------------------------------
mtcars2$transmission <-
  factor(mtcars2$am, levels = c(0, 1), labels = c("Automatic", "Manual"))

## -----------------------------------------------------------------------------
mtcars2$make  <- sub("^(\\w+)\\s(.+)", "\\1", rownames(mtcars2))
mtcars2$model <- sub("^(\\w+)\\s(.+)", "\\2", rownames(mtcars2))
rownames(mtcars2) <- NULL

## -----------------------------------------------------------------------------
set.seed(42)
mtcars2$test_date <-
  as.POSIXct("1974-01-03", tz = "GMT") +
  cumsum(sample(c(2, 3, 4, 7) * 3600 * 24, size = nrow(mtcars2), replace = TRUE))

## -----------------------------------------------------------------------------
mtcars2 <-
  mtcars2[, c("make", "model", "mpg", "disp", "hp", "drat", "wt", "qsec",
              "cyl", "cyl_character", "cyl_factor",
              "vs", "engine",
              "am", "transmission",
              "gear", "gear_factor",
              "carb",
              "test_date")]

## -----------------------------------------------------------------------------
str(mtcars2)

## -----------------------------------------------------------------------------
# copied text from the manuscript
pefr_table <-
  read.delim(
             header = FALSE,
             text = "
1	494	490	512	525
2	395	397	430	415
3	516	512	520	508
4	434	401	428	444
5	476	470	500	500
6	557	611	600	625
7	413	415	364	460
8	442	431	380	390
9	650	638	658	642
10	433	429	445	432
11	417	420	432	420
12	656	633	626	605
13	267	275	260	227
14	478	492	477	467
15	178	165	259	268
16	423	372	350	370
17	427	421	451	443")

## -----------------------------------------------------------------------------
pefr <-
  expand.grid(subject = 1:17,
              measurement = 1:2,
              meter   = c("Wright peak flow meter", "Mini Wright peak flow meter"),
              KEEP.OUT.ATTRS = FALSE,
              stringsAsFactors = FALSE)
pefr$pefr <- do.call(c, pefr_table[, 2:5])

head(pefr)

## ----label = "sessioninfo"----------------------------------------------------
sessionInfo()

