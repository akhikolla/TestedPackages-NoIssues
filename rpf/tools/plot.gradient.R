library(rpf)
library(ggplot2)
library(reshape2)
library(testthat)
#options(digits=22)
#options(width=200)

# need to integrate TODO
inspect.gradient <- function(spec, param) {
  w <- seq(1, spec@numOutcomes) * 5.0
  samples <- 50
  got <- list()
  for (px in 1:length(param)) {
    orig <- param[px]
    df1 <- data.frame(grid=orig + seq(-.5, .5, length.out=samples+1))
    for (gx in 1:length(df1$grid)) {
      param[px] <- df1$grid[gx]
      df1$LL[gx] <- sum(w * rpf.logprob(spec, param, loc))
    }
    df2 <- data.frame(empirical=diff(df1$LL) / diff(df1$grid),
                      grid=df1$grid[1:samples] + diff(df1$grid)/2)
    for (gx in 1:length(df2$grid)) {
      param[px] <- df2$grid[gx]
      df2$analytic[gx] <- rpf.gradient(spec, param, loc, w)[px]
    }
    df2$px <- px
    param[px] <- orig
    got <- rbind(got, df2)
  }
  got$px <- factor(got$px)
  as.data.frame(got)
}

plot.gradient <- function(df, px) {
  df2 <- df[df$px == px,]
  df2 <- melt(df2, id.vars=c('px', 'grid'), value.name="dLL", variable.name="type")
  spot <- df2[quantile(1:dim(df2)[[1]], probs=seq(0,1,.2)),]
  
  ggplot(df2, aes(grid, dLL, color=type)) + geom_line() +
    geom_point(data=spot, aes(x=grid,y=dLL), alpha=.5, size=3, color='yellow')
}

plot.gradients <- function(df, px) {
  df2 <- df[df$px == px,]
  df2 <- df2[,c('empirical', 'grid')]
  colnames(df2)[[1]] <- "dLL"
  df2$name <- paste("empirical", px, sep="")
  df3 <- df[,c('analytic', 'grid', 'px')]
  colnames(df3)[[1]] <- "dLL"
  df3$name <- paste("analytic", df3$px, sep="")
  df3$px <- NULL
  plot.data <- rbind(df2, df3)
  plot.data$name <- factor(plot.data$name)
  spot <- df2[quantile(1:dim(df2)[[1]], probs=seq(0,1,.2)),]
  ggplot(plot.data, aes(grid, dLL, color=name)) + geom_line() +
    geom_point(data=spot, aes(x=grid,y=dLL), alpha=.5, size=3, color='yellow')
}

spec <- list()
param <- list()
spec [[length(spec) +1]] <- rpf.drm(poor=TRUE)
param[[length(param)+1]] <- c(1, 0, 0)
spec [[length(spec) +1]] <- rpf.drm(multidimensional=TRUE)
param[[length(param)+1]] <- c(1, 0, 0)
spec [[length(spec) +1]] <- rpf.gpcm(3)
param[[length(param)+1]] <- c(1, -10, 10)

#for (mx in 1:length(spec)) {
  mx<-3
  ispec <- spec[[mx]]
  iparam <- param[[mx]]
  df <- inspect.gradient(ispec, iparam)
  expect_equal(df$empirical, df$analytic, tolerance = 1e-4)
#}

# diagnostic tools
if (0) {
  plot.gradient(df, 2)
  plot.gradients(df, 1)
}
