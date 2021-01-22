set.seed(2)
library(data.table)
signal <- c(
  rnorm(25, mean = 10),
  rnorm(25, mean = 7),
  rnorm(25, mean = 8),
  rnorm(25, mean = 5))
#outliers
signal[86] <- 10
labels.dt <- data.table(
  start = c(20, 45, 80),
  end = c(30, 55, 90),
  changes = c(1, 1, 0))
signal.dt <- data.table(
  signal,
  position=seq_along(signal))
label.colors <- c(
  "1"="#ff7d7d",
  "0"="#f6c48f")
sig.color <- "grey50"
if(require(ggplot2)){
  gg.data <- ggplot()+
    geom_rect(aes(
      xmin=start, xmax=end,
      fill=paste(changes),
      ymin=-Inf, ymax=Inf),
      alpha=0.5,
      data=labels.dt)+
    geom_point(aes(
      position, signal),
      color=sig.color,
      data=signal.dt)+
    scale_x_continuous(
      "position",
      breaks=seq(0, 100, by=10))+
    scale_fill_manual("label", values=label.colors)+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))
  print(gg.data)
}

label.list <- list(
  OPART=labels.dt[0],
  LOPART=labels.dt)
seg.dt.list <- list()
change.dt.list <- list()
cost.dt.list <- list()
for(model.name in names(label.list)){
  label.dt <- data.table(label.list[[model.name]])
  fit <- LOPART::LOPART(signal, label.dt, 10)
  Algorithm <- factor(model.name, names(label.list))
  tau.dt <- fit$cost[, .(
    cost_candidates,
    tau=0:(.N-1),
    change=seq_along(cost_candidates)-0.5
  )]
  cost.dt.list[[model.name]] <- data.table(Algorithm, tau.dt)
  seg.dt.list[[model.name]] <- data.table(Algorithm, fit$segments)
  change.dt.list[[model.name]] <- data.table(Algorithm, fit$changes)
}
seg.dt <- do.call(rbind, seg.dt.list)
change.dt <- do.call(rbind, change.dt.list)
cost.dt <- do.call(rbind, cost.dt.list)

algo.sizes <- c(
  OPART=1,
  LOPART=0.5)
algo.colors <- c(
  OPART="deepskyblue",
  LOPART="black")
algo.shapes <- c(
  OPART=1,
  LOPART=2)
if(require(ggplot2)){
  gg.data+
    scale_size_manual(values=algo.sizes)+
    scale_color_manual(values=algo.colors)+      
    geom_vline(aes(
      xintercept=change,
      size=Algorithm,
      color=Algorithm),
      data=change.dt)+
    geom_segment(aes(
      start-0.5, mean,
      size=Algorithm,
      color=Algorithm,
      xend=end+0.5, yend=mean),
      data=seg.dt)
}

if(require(ggplot2)){
  ggplot()+
    geom_rect(aes(
      xmin=start, xmax=end,
      fill=paste(changes),
      ymin=-Inf, ymax=Inf),
      alpha=0.5,
      data=labels.dt)+
    scale_fill_manual("label", values=label.colors)+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    scale_x_continuous(
      "position",
      breaks=seq(0, 100, by=10))+
    geom_point(aes(
      change, cost_candidates,
      color=Algorithm, shape=Algorithm),
      data=cost.dt)+
    scale_color_manual(values=algo.colors)+      
    scale_shape_manual(values=algo.shapes)
}

abbrev.vec <- c(
  data="data and models",
  cost="cost of last change")
yfac <- function(l){
  factor(abbrev.vec[[l]], abbrev.vec)
}
COST <- function(dt){
  data.table(y.var=yfac("cost"), dt)
}
DATA <- function(dt){
  data.table(y.var=yfac("data"), dt)
}
if(require(ggplot2)){
  ggplot()+
    geom_rect(aes(
      xmin=start, xmax=end,
      fill=paste(changes),
      ymin=-Inf, ymax=Inf),
      alpha=0.5,
      data=labels.dt)+
    scale_fill_manual("label", values=label.colors)+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(y.var ~ ., scales="free")+
    geom_vline(aes(
      xintercept=change,
      size=Algorithm,
      color=Algorithm),
      data=change.dt)+
    geom_segment(aes(
      start-0.5, mean,
      size=Algorithm,
      color=Algorithm,
      xend=end+0.5, yend=mean),
      data=DATA(seg.dt))+
    geom_point(aes(
      position, signal),
      color=sig.color,
      shape=1,
      data=DATA(signal.dt))+
    scale_size_manual(values=algo.sizes)+
    scale_color_manual(values=algo.colors)+
    scale_shape_manual(values=algo.shapes)+
    ylab("")+
    scale_x_continuous(
      "position",
      breaks=seq(0, 100, by=10))+
    geom_point(aes(
      change, cost_candidates,
      color=Algorithm, shape=Algorithm),
      data=COST(cost.dt))
}

