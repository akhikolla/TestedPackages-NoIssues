---
output:
  html_document: default
  pdf_document: default
---
# YPInterimTesting

Provide monitoring boundaries and nominal p-values at the interim looks using the adaptively weighted log-rank test developed by Yang and Prentice (2010). The package use a re-sampling method to obtain stopping boundaries in sequential designs. The asymptotic distribution of the test statistics of the adaptively weighted log-rank test at the interim looks is examined in Yang (2018).

## Installation
``` r
install.packages("YPInterimTesting")
```


## Example

```
library(YPInterimTesting)
data(virtual_data) # the data created to show how to utilize the package

time <- virtual$time
event <- virtual$event
group <- virtual$group

spendfun <- c(1.3E-5, 4.4E-4, 0.003, 0.008)

result_all <- ypinterim(time, event, group, spendfun=spendfun)
result_all
summary(result_all)

```
The above example shows how to test the package with a historical data where interim data at all looks are available. The object `result_all` can be formatted to a table using the function `summary`. 

```
summary(result)
```
The known boundaries of the first few looks can be supplied using `crtivalue`. The following three examples show how to supply `crtivalue` when the boundaries at the previous looks exist.

We first need to calculate the boundary at the first look. The spending function value at the first look is needed.

```
# Assume that we only have the information on the first look
time <- virtual$time[, 1]
event <- virtual$event[, 1]
group <- virtual$group

spendfun <- c(1.3E-5)

result_look1 <- ypinterim(time, event, group, spendfun=spendfun)
result_look1
summary(result_look1)
```

When calculating the boundary at the second look, the spending function at the two looks, and boundary at the first look (i.e., the value obtained from the previous example), should be supplied.

```
time <- virtual$time[, 1:2]
event <- virtual$event[, 1:2]
group <- virtual$group

spendfun <- c(1.3E-5, 4.4E-4)
critvalue <- c(4.36) # the boundary of the first look is supplied.

result_look2 <- ypinterim(time, event, group, spendfun=spendfun, critvalue = critvalue)
result_look2
summary(result_look2)
```

Similarly, when calculating the boundary at the third look, the spending function at the three looks, and boundaries at the first two looks, should be supplied.

```
time <- virtual$time[, 1:3]
event <- virtual$event[, 1:3]
group <- virtual$group

spendfun <- c(1.3E-5, 4.4E-4, 0.003)
critvalue <- c(4.36, 3.42) # the boundaries at the first two looks are supplied.

result_look3 <- ypinterim(time, event, group, spendfun=spendfun, critvalue = critvalue)
result_look3
summary(result_look3)
```

## Reference
Yang, S. (2018). Interim monitoring using the adaptively weighted log-rank test in clinical trials for survival outcomes. Statistics in Medicine.

Yang, S., & Prentice, R. (2010). Improved logrank-type tests for survival data using adaptive weights. Biometrics, 66(1), 30-38.

