## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE-------------------------------------------------------
#  # install package (to be done once)
#  install.packages(EloChoice)

## ------------------------------------------------------------------------
# load package (every time you want to use the package)
library(EloChoice) 

## ---- eval=FALSE---------------------------------------------------------
#  # Windows
#  xdata <- read.table(file = "c:\\datafiles\\myfile.txt", sep = "\t", header = TRUE)
#  
#  # Mac
#  xdata <- read.table(file = "/Volumes/mydrive/myfile.txt", sep = "\t", header = TRUE)
#  
#  str(xdata)

## ----exampletab, echo=FALSE, results='markdown'--------------------------
# invisible(table_nums(name = "exampletab", caption = "A"))

winner <- c("ab", "cf", "ab", "dd", "ab")
loser <- c("cf", "xs", "xs", "cf", "cf")
rater <- c("A", "A", "A", "A", "B")
date <- c("2010-01-01", "2010-01-01", "2010-01-01", "2010-01-01", "2010-01-04")
time <- c("14:34:01", "14:34:08", "14:34:11", "14:34:15", "09:17:20")
mytab <- data.frame(winner, loser, rater, date, time)
colnames(mytab) <- c("preferred stimulus", "losing stimulus", "rater", "date", "time")
cap <- "A possible data set layout. Note that R replaces spaces in column names with periods during the reading step."
knitr::kable(mytab, caption = cap)

## ------------------------------------------------------------------------
set.seed(123)
xdata <- randompairs(nstim = 7, nint = 700, reverse = 0.1)
head(xdata)

## ------------------------------------------------------------------------
set.seed(123)
res <- elochoice(winner = xdata$winner, loser = xdata$loser, runs = 1000)
summary(res)

## ------------------------------------------------------------------------
ratings(res, show = "original", drawplot = FALSE)

## ------------------------------------------------------------------------
ratings(res, show = "mean", drawplot = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  myratings <- ratings(res, show = "mean", drawplot = FALSE)
#  # Windows
#  xdata <- write.table(myratings, "c:\\datafiles\\myratings.txt", sep = "\t", header = TRUE)
#  # Mac
#  xdata <- write.table(myratings, "/Volumes/mydrive/myratings.txt", sep = "\t", header = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  myratings <- ratings(res, show = "all", drawplot = FALSE)
#  # Windows
#  xdata <- write.table(myratings, "c:\\datafiles\\myratings.txt", sep = "\t", header = TRUE, row.names = FALSE)
#  # Mac
#  xdata <- write.table(myratings, "/Volumes/mydrive/myratings.txt", sep = "\t", header = TRUE, row.names = FALSE)

## ----fig1plot, fig.align='center', fig.cap="Elo ratings of 7 stimuli after 700 rating events and 1000 randomizations of the sequence. The black circles represent the average rating at the end of the 1000 generated sequences for each stimulus, and the black lines represent their ranges. The grey circles show the final ratings from the original sequence.", echo = 2:2, fig.width=7, fig.height=4----
par(mar = c(4.1, 4.1, 0.5, 0.5), family = "serif")
ratings(res, show = NULL, drawplot = TRUE)

## ----upset, echo=FALSE, results='markdown'-------------------------------
pref <- c(1, 1, 0, 0, 1, 0, 1, 0, 0, 0)
upset <- c("yes", "yes", "no", "no", "yes", "no", "yes", "no", "no", "no" )
ratingdiff  <- c(200, 300, 100, 150, 200, 140, 280, 90, 150, 120)
ratingdiff2 <- c(90, 100, 300, 280, 120, 200, 140, 150, 200, 150)
mytab <- data.frame(pref, upset, ratingdiff, ratingdiff2)
colnames(mytab) <- c("higher rated is not preferred", "upset", "rating difference (example 1)", "rating difference (example 2)")

cap <- "10 rating decisions that were either in accordance with the expectation or not. Two different rating differences are given to illustrate the weighted upset index. Note that the values are the same, just their assignment to different interactions is changed and consequently the column means are the same for both."
knitr::kable(mytab, caption = cap)

## ------------------------------------------------------------------------
upsets <- reliability(res)
head(upsets)

## ------------------------------------------------------------------------
mean(upsets$upset)
mean(upsets$upset.wgt)

## ------------------------------------------------------------------------
set.seed(123)
xdata <- randompairs(nstim = 7, nint = 700, reverse = 0.3)
res <- elochoice(winner = xdata$winner, loser = xdata$loser, runs = 1000)
upsets <- reliability(res)
mean(upsets$upset)
mean(upsets$upset.wgt)

## ------------------------------------------------------------------------
data(physical)
# limit to 10 raters
physical <- subset(physical, raterID %in% c(1, 2, 8, 10, 11, 12, 23, 27, 31, 47))
set.seed(123)
res <- raterprog(physical$Winner, physical$Loser, physical$raterID, progbar = FALSE)

## ----raterprog1, fig.width=7, fig.height = 4, fig.align='center', fig.cap="Reliability index $R'$ as function of number of raters included in the rating process. In this example, it appears as if including 5 raters is sufficient, as additional raters improve rating reliability relatively little.", echo = 2:2----
par(mar = c(4.1, 4.1, 0.5, 0.5), family = "serif")
raterprogplot(res)

## ------------------------------------------------------------------------
set.seed(123)
res <- raterprog(physical$Winner, physical$Loser, physical$raterID, progbar = FALSE, ratershuffle = 10)

## ----raterprog2, fig.width=7, fig.height = 4, fig.align='center', fig.cap="Reliability index $R'$ as function of number of raters included in the rating process. Here, we used the original rater order and an additional 9 random orders. Grey bars reflect quartiles, grey points are the results from the original rater order (compare to figure~\ref{fig:raterprog1}), and black points are the average values of the 10 rater orders used. The interpretation would likely to be the same as for figure~\ref{fig:raterprog1}.", echo=2:2----
par(mar = c(4.1, 4.1, 0.5, 0.5), family = "serif")
raterprogplot(res)

## ------------------------------------------------------------------------
data(physical)
set.seed(123)
res <- elochoice(winner = physical$Winner, loser = physical$Loser, runs = 500)
summary(res)
ratings(res, show = "mean", drawplot = FALSE)

## ----fig2plot, fig.width=7, fig.height = 4, fig.align='center', fig.cap="Elo ratings of 82 stimuli after 4,592 rating events and 500 randomizations of the sequence. The black circles represent the average (mean) rating at the end of the 500 generated sequences, and the black lines represent their ranges. The grey circles show the final ratings from the original sequence. Note that not all stimulus IDs fit on the x-axis, so most are omitted.", echo = 2:2----
par(mar = c(4.1, 4.1, 0.5, 0.5), family = "serif")
ratings(res, show = NULL, drawplot = TRUE)

## ------------------------------------------------------------------------
# total of seven trials with two 'self-trials' (trials 6 and 7)
w <- c(letters[1:5], "a", "b"); l <- c(letters[2:6], "a", "b")
res <- elochoice(w, l)
ratings(res, drawplot=FALSE)
summary(res)
# total of five trials without 'self-trials'
w <- c(letters[1:5]); l <- c(letters[2:6])
res <- elochoice(w, l)
ratings(res, drawplot=FALSE)
summary(res)

