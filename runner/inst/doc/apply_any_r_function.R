## ----eval=FALSE---------------------------------------------------------------
#  # full windows
#  runner(1:15)
#  
#  # summarizing - sum
#  runner(
#    1:15,
#    f = sum
#  )
#  
#  # summarizing - concatenating
#  runner(
#    letters[1:15],
#    f = paste,
#    collapse = " > "
#  )

## ----eval=FALSE---------------------------------------------------------------
#  # summarizing - sum of 4-elements
#  runner(
#    1:15,
#    k = 4,
#    f = sum
#  )
#  
#  # summarizing - slope from lm
#  df <- data.frame(
#    a = 1:15,
#    b = 3 * 1:15 + rnorm(15)
#  )
#  
#  runner(
#    x = df,
#    k = 5,
#    f = function(x) {
#      model <- lm(b ~ a, data = x)
#      coefficients(model)["a"]
#    }
#  )

## ----eval=FALSE---------------------------------------------------------------
#  idx <- c(4, 6, 7, 13, 17, 18, 18, 21, 27, 31, 37, 42, 44, 47, 48)
#  
#  # summarize - mean
#  runner::runner(
#    x = idx,
#    k = 5, # 5-days window
#    lag = 1,
#    idx = idx,
#    f = function(x) mean(x)
#  )
#  
#  
#  # use Date or datetime sequences
#  runner::runner(
#    x = idx,
#    k = "5 days", # 5-days window
#    lag = 1,
#    idx = Sys.Date() + idx,
#    f = function(x) mean(x)
#  )
#  
#  # obtain window from above illustration
#  runner::runner(
#    x = idx,
#    k = "5 days",
#    lag = 1,
#    idx = Sys.Date() + idx
#  )

## ----eval=FALSE---------------------------------------------------------------
#  idx <- c(4, 6, 7, 13, 17, 18, 18, 21, 27, 31, 37, 42, 44, 47, 48)
#  
#  # summary
#  runner::runner(
#    x = 1:15,
#    k = 5,
#    lag = 1,
#    idx = idx,
#    at = c(18, 27, 48, 31),
#    f = mean
#  )
#  
#  # full window
#  runner::runner(
#    x = idx,
#    k = 5,
#    lag = 1,
#    idx = idx,
#    at = c(18, 27, 48, 31)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  idx_date <- seq(Sys.Date(), Sys.Date() + 365, by = "1 month")
#  
#  # change interval to 4-months
#  runner(
#    x = 0:12,
#    idx = idx_date,
#    at = "4 months"
#  )
#  
#  # calculate correlation at every 6-months
#  runner(
#    x = data.frame(
#      a = 1:13,
#      b = 1:13 + rnorm(13, sd = 5),
#      idx_date
#    ),
#    idx = "idx_date",
#    at = "6 months",
#    f = function(x) {
#      cor(x$a, x$b)
#    }
#  )

## ----eval=FALSE---------------------------------------------------------------
#  # summarizing - concatenating
#  runner::runner(
#    x = 1:10,
#    lag = c(-1, 2, -1, -2, 0, 0, 5, -5, -2, -3),
#    k = c(0, 1, 1, 1, 1, 5, 5, 5, 5, 5),
#    f = paste,
#    collapse = ","
#  )
#  
#  # full window
#  runner::runner(
#    x = 1:10,
#    lag = 1,
#    k = c(1, 1, 1, 1, 1, 5, 5, 5, 5, 5)
#  )
#  
#  # on dates
#  idx <- c(4, 6, 7, 13, 17, 18, 18, 21, 27, 31, 37, 42, 44, 47, 48)
#  
#  runner::runner(
#    x = 1:15,
#    lag = sample(c("-2 days", "-1 days", "1 days", "2 days"),
#                 size = 15,
#                 replace = TRUE),
#    k = sample(c("5 days", "10 days", "15 days"),
#               size = 15,
#               replace = TRUE),
#    idx = Sys.Date() + idx,
#    f = function(x) mean(x)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  idx <- c(4, 6, 7, 13, 17, 18, 18, 21, 27, 31, 37, 42, 44, 47, 48)
#  
#  runner::runner(
#    x = 1:15,
#    k = 5,
#    lag = 1,
#    idx = idx,
#    at = c(4, 18, 48, 51),
#    na_pad = TRUE,
#    f = function(x) mean(x)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  x <- cumsum(rnorm(40))
#  y <- 3 * x + rnorm(40)
#  date <- Sys.Date() + cumsum(sample(1:3, 40, replace = TRUE)) # unequaly spaced time series
#  group <-  rep(c("a", "b"), 20)
#  
#  df <- data.frame(date, group, y, x)
#  
#  slope <- runner(
#    df,
#    function(x) {
#      coefficients(lm(y ~ x, data = x))[2]
#    }
#  )
#  
#  plot(slope)

## ----eval=FALSE---------------------------------------------------------------
#  library(dplyr)
#  
#  df <- df %>%
#    group_by(group) %>%
#    mutate(
#      cumulative_mse = runner(
#        x = .,
#        k = "20 days",
#        idx = "date", # specify column name instead df$date
#        f = function(x) {
#          coefficients(lm(y ~ x, data = x))[2]
#        }
#      )
#    )
#  
#  library(ggplot2)
#  df %>%
#    ggplot(aes(x = date, y = cumulative_mse, group = group, color = group)) +
#    geom_line()
#  

## ----eval=FALSE---------------------------------------------------------------
#  df %>%
#    group_by(group) %>%
#    run_by(idx = "date", k = "20 days", na_pad = FALSE) %>%
#    mutate(
#      cumulative_mse = runner(
#        x = .,
#        f = function(x) {
#          mean((residuals(lm(y ~ x, data = x))) ^ 2)
#        }
#      ),
#  
#      intercept = runner(
#        x = .,
#        f = function(x) {
#          coefficients(lm(y ~ x, data = x))[1]
#        }
#      ),
#  
#      slope = runner(
#        x = .,
#        f = function(x) {
#          coefficients(lm(y ~ x, data = x))[2]
#        }
#      )
#    )
#  
#  

