## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(crandep)
library(dplyr)
library(igraph)

## -----------------------------------------------------------------------------
get_dep("dplyr", "Imports")
get_dep("MASS", "depends")

## -----------------------------------------------------------------------------
get_dep_df("dplyr", c("imports", "LinkingTo"))

## -----------------------------------------------------------------------------
get_dep("abc", "depends")
get_dep("abc", "reverse_depends")
get_dep_df("abc", c("depends", "reverse_depends"))

## ---- echo=FALSE--------------------------------------------------------------
data.frame(from = "A", to = "B", type = "c", reverse = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
data.frame(from = "B", to = "A", type = "c", reverse = TRUE)

## -----------------------------------------------------------------------------
df0.abc <- get_dep_df("abc", "all")
df0.abc
df0.rstan <- get_dep_df("rstan", "all")
dplyr::count(df0.rstan, type, reverse) # all 8 types

## ----echo=FALSE---------------------------------------------------------------
df0.all <- get_dep_all_packages()
v0.all <- df0.all %>% group_by(from) %>% count(a = n_distinct(type, reverse)) %>% filter(a == 8L)

## -----------------------------------------------------------------------------
df0.imports <- rbind(
    get_dep_df("ggplot2", "Imports"),
    get_dep_df("dplyr", "Imports"),
    get_dep_df("tidyr", "Imports"),
    get_dep_df("readr", "Imports"),
    get_dep_df("purrr", "Imports"),
    get_dep_df("tibble", "Imports"),
    get_dep_df("stringr", "Imports"),
    get_dep_df("forcats", "Imports")
)
head(df0.imports)
tail(df0.imports)

## ---- out.width="660px", out.height="660px", fig.width=12, fig.height=12, fig.show="hold"----
g0.imports <- igraph::graph_from_data_frame(df0.imports)
set.seed(1457L)
old.par <- par(mar = rep(0.0, 4))
plot(g0.imports, vertex.label.cex = 1.5)
par(old.par)

## -----------------------------------------------------------------------------
igraph::is_dag(g0.imports)

## ---- out.width="660px", out.height="660px", fig.width=12, fig.height=12, fig.show="hold"----
df0.nodes <- data.frame(name = c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats"), stringsAsFactors = FALSE)
g0.core <- df_to_graph(df0.imports, df0.nodes)
set.seed(259L)
old.par <- par(mar = rep(0.0, 4))
plot(g0.core, vertex.label.cex = 1.5)
par(old.par)

## -----------------------------------------------------------------------------
topo_sort_kahn(g0.core)

## -----------------------------------------------------------------------------
set.seed(387L); topo_sort_kahn(g0.core, random = TRUE)

## -----------------------------------------------------------------------------
df0.topo <- topo_sort_kahn(g0.imports)
head(df0.topo)
tail(df0.topo)

