## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(crandep)
library(dplyr)
library(ggplot2)
library(igraph)
library(visNetwork)

## -----------------------------------------------------------------------------
data(cran_dependencies)
cran_dependencies
dplyr::count(cran_dependencies, type, reverse)

## -----------------------------------------------------------------------------
df0.cran <- get_dep_all_packages()
head(df0.cran)
dplyr::count(df0.cran, type, reverse) # numbers in general larger than above

## -----------------------------------------------------------------------------
g0.depends <- get_graph_all_packages(type = "depends")
g0.rev_depends <- get_graph_all_packages(type = "reverse depends")
g0.depends
g0.rev_depends

## -----------------------------------------------------------------------------
g1.depends <- df0.cran %>%
    dplyr::filter(type == "depends" & !reverse) %>%
    df_to_graph(nodelist = dplyr::rename(df0.cran, name = from))
g1.rev_depends <- df0.cran %>%
    dplyr::filter(type == "depends" & reverse) %>%
    df_to_graph(nodelist = dplyr::rename(df0.cran, name = from))
g1.depends # same as g0.depends
g1.rev_depends # same as g0.rev_depends

## -----------------------------------------------------------------------------
igraph::is_dag(g0.depends)
igraph::is_dag(g0.rev_depends)

## -----------------------------------------------------------------------------
df1.rev_depends <- df0.cran %>%
    dplyr::filter(type == "depends" & reverse) %>%
    df_to_graph(nodelist = NULL, gc = FALSE) %>%
    igraph::as_data_frame() # to obtain the edge list
df1.depends <- df0.cran %>%
    dplyr::filter(type == "depends" & !reverse) %>%
    df_to_graph(nodelist = NULL, gc = FALSE) %>%
    igraph::as_data_frame()
dfa.diff.depends <- dplyr::anti_join(
    df1.rev_depends,
    df1.depends,
    c("from" = "to", "to" = "from")
)
head(dfa.diff.depends)

## -----------------------------------------------------------------------------
dfb.diff.depends <- dplyr::anti_join(
    df1.depends,
    df1.rev_depends,
    c("from" = "to", "to" = "from")
)
head(dfb.diff.depends)

## -----------------------------------------------------------------------------
df0.summary <- dplyr::count(df0.cran, from, type, reverse)
head(df0.summary)

## -----------------------------------------------------------------------------
df0.summary %>%
    dplyr::filter(reverse) %>%
    dplyr::group_by(type) %>%
    dplyr::top_n(1, n)

## ---- out.width="660px", out.height="660px", fig.width=9, fig.height=9--------
df1.summary <- df0.summary %>%
    dplyr::count(type, reverse, n)
gg0.summary <- df1.summary %>%
    dplyr::mutate(reverse = ifelse(reverse, "reverse", "forward")) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(n, nn)) +
    ggplot2::facet_grid(type ~ reverse) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Degree", y = "Number of packages") +
    ggplot2::theme_bw(20)
gg0.summary

## -----------------------------------------------------------------------------
prefix <- "http://CRAN.R-project.org/package=" # canonical form
degrees <- igraph::degree(g0.depends)
set.seed(3579L)
clusters <- igraph::membership(igraph::cluster_spinglass(g0.depends))
df0.nodes <- dplyr::full_join(
    data.frame(id = names(degrees), value = degrees),
	data.frame(id = names(clusters), group = as.character(clusters)),
	"id") %>%
    dplyr::mutate(title = paste0('<a href=\"', prefix, id, '\">', id, '</a>'))
df0.edges <- igraph::as_data_frame(g0.depends, what = "edges")

## -----------------------------------------------------------------------------
set.seed(2345L)
vis0 <- visNetwork::visNetwork(df0.nodes, df0.edges, width = "100%", height = "720px") %>%
    visNetwork::visOptions(highlightNearest = TRUE) %>%
	visNetwork::visEdges(arrows = "to", color = list(opacity = 0.5)) %>%
	visNetwork::visNodes(fixed = TRUE) %>%
	visNetwork::visIgraphLayout(layout = "layout_with_drl")
vis0

