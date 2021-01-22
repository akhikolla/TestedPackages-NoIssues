## ----label = "setup", include = FALSE-----------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

## -----------------------------------------------------------------------------
library(qwraps2)
packageVersion("qwraps2")
library(ggplot2)

## -----------------------------------------------------------------------------
set.seed(42)
n <- 250
x1 <- x2 <- x3 <- x4 <- vector('numeric', length = n)
x1[1] <- runif(1)
x2[1] <- runif(1)
x3[1] <- runif(1)
x4[1] <- runif(1)

# white noise
Z.1 <- rnorm(n, 0, 1)
Z.2 <- rnorm(n, 0, 2)
Z.3 <- rnorm(n, 0, 5)

for(i in 2:n)
{
  x1[i] <- x1[i-1] + Z.1[i] - Z.1[i-1] + x4[i-1] - x2[i-1]
  x2[i] <- x2[i-1] - 2 * Z.2[i] + Z.2[i-1] - x4[i-1]
  x3[i] <- x3[i-1] + x2[i-1] + 0.2 * Z.3[i] + Z.3[i-1]
  x4[i] <- x4[i-1] + runif(1, 0.5, 1.5) * x4[i-1]
}
testdf <- data.frame(x1, x2, x3, x4)

# Base acf plot for one variable
acf(testdf$x1)

# qacf plot for one variable
qacf(testdf$x1)
qacf(testdf$x1, show_sig = TRUE)


## ----fig.width = 5, fig.height = 5--------------------------------------------
# more than one variable
acf(testdf)
qacf(testdf)
qacf(testdf, show_sig = TRUE)

## -----------------------------------------------------------------------------
acf_plot_data <- qacf(testdf)$data
head(acf_plot_data)

## -----------------------------------------------------------------------------
pefr_m1 <-
  cbind("Large" = pefr[pefr$measurement == 1 & pefr$meter == "Wright peak flow meter", "pefr"],
        "Mini"  = pefr[pefr$measurement == 1 & pefr$meter == "Mini Wright peak flow meter", "pefr"])

## -----------------------------------------------------------------------------
cor(pefr_m1)

qplot(x = pefr_m1[, 1],
      y = pefr_m1[, 2],
      geom = "point",
      xlab = "Large Meter",
      ylab = "Mini Meter",
      xlim = c(0, 800),
      ylim = c(0, 800)) +
geom_abline(slope = 1)

## -----------------------------------------------------------------------------
qblandaltman(pefr_m1) +
xlim(0, 800) +
ylim(-100, 100) +
xlab("Average of two meters") +
ylab("Difference in the measurements")

## -----------------------------------------------------------------------------
pefr_mini <-
  cbind(m1 = pefr[pefr$measurement == 1 & pefr$meter == "Mini Wright peak flow meter", "pefr"],
        m2 = pefr[pefr$measurement == 2 & pefr$meter == "Mini Wright peak flow meter", "pefr"])

qblandaltman(pefr_mini)

## -----------------------------------------------------------------------------
# create a survfit object
require(survival)
leukemia.surv <- survival::survfit(survival::Surv(time, status) ~ x, data = survival::aml)

# base R km plot
survival:::plot.survfit(leukemia.surv, conf.int = TRUE, lty = 2:3, col = 1:2)


## ----fig.width = 5------------------------------------------------------------
# qkmplot
qkmplot(leukemia.surv, conf_int = TRUE)

# build a data.frame for plotting km curves, this could be helpful for
# creating bespoke plots
leukemia_km_data <- qkmplot_bulid_data_frame(leukemia.surv)
head(leukemia_km_data, 3)


## ----fig.width = 5------------------------------------------------------------
qkmplot(leukemia_km_data)


## ----fig.width = 5------------------------------------------------------------
# intercept only plot
intonly_fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = survival::aml)
survival:::plot.survfit(intonly_fit, conf.int = TRUE)
qkmplot(intonly_fit, conf_int = TRUE)

## -----------------------------------------------------------------------------
data(diamonds, package = "ggplot2")

# Create two logistic regression models
fit1 <- glm(I(price > 2800) ~ cut * color, data = diamonds, family = binomial())
fit2 <- glm(I(price > 2800) ~ cut + color + clarity, data = diamonds, family = binomial())

# Easiest way to get an ROC plot:
qroc(fit1)
qroc(fit2)

# Create two data sets, this will also let you get the AUC out
data1 <- qroc_build_data_frame(fit1)
data2 <- qroc_build_data_frame(fit2)

auc(data1)
auc(data2)

# Plotting the ROC from the data set can be done too
qroc(data1)

# Add the AUC value to the plot title
qroc(data2) + ggtitle(paste("Fit 2\nAUC =", round(auc(data2), 2)))

# build a data set for plotting to ROCs on one plot
plot_data <- rbind(cbind(Model = "fit1", data1),
                   cbind(Model = "fit2", data2))
qroc(plot_data) + aes(color = Model)

# with AUC in the legend
plot_data <- rbind(cbind(Model = paste("Fit1\nauc =", round(auc(data1), 3)), data1),
                   cbind(Model = paste("Fit2\nauc =", round(auc(data2), 3)), data2))
qroc(plot_data) +
  theme_bw() +
  aes(color = Model, linetype = Model) +
  theme(legend.position   = "bottom",
        legend.text.align = 0.5)

## ----label = "sessioninfo"----------------------------------------------------
sessionInfo()

