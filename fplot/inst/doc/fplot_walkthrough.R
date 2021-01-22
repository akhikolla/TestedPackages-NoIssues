## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, include = FALSE-------------------------------------------
# We conditionnaly run the code pieces
IS_FIXEST = TRUE
msg = function() NULL
if(!requireNamespace("fixest", quietly = TRUE)){
  IS_FIXEST = FALSE
	message("This vignette uses a data set from the package 'fixest', which is not currently available. Examples relating to this data set are removed.")
	msg = function() message("Package 'fixest' not available: result not reported.")
}
library(fplot)

## -----------------------------------------------------------------------------
setFplot_dict(c(Origin = "Exporting Country", Destination = "Importing Country", Euros = "Exports Value in €", jnl_top_25p = "Pub. in Top 25% journal", jnl_top_5p = "Publications in Top 5% journal", journal = "Journal", institution = "U.S. Institution", Petal.Length = "Petal Length"))

## ---- echo = TRUE-------------------------------------------------------------
library(fplot)
data(us_pub_econ, package = "fplot")

## ---- echo=FALSE, results='asis'----------------------------------------------
tab = head(us_pub_econ)
# Following line to be run when pdf output only
if ("pdf_document" %in% rmarkdown::all_output_formats(knitr::current_input())){
  tab = as.data.table(lapply(tab, fplot:::truncate_string, method = "mid", trunc = 16))
}
knitr::kable(tab)

## ---- fig.width=7-------------------------------------------------------------
# When there is only one variable, you can use a vector
plot_distr(us_pub_econ$institution)


## ---- fig.width=7-------------------------------------------------------------
plot_distr(jnl_top_5p ~ institution, us_pub_econ)


## ---- fig.width=7-------------------------------------------------------------
plot_distr(1 + jnl_top_25p + jnl_top_5p ~ institution, us_pub_econ)


## ---- fig.width=7-------------------------------------------------------------
plot_distr(~ journal | institution, us_pub_econ)


## ---- fig.width=7-------------------------------------------------------------
# Previous graph + "other" column
plot_distr(~ journal | institution, us_pub_econ, other = TRUE)


## ---- fig.width=7-------------------------------------------------------------
plot_distr(~ journal | institution, us_pub_econ, mod.select = "harvard|boston.+uni|mass.+inst")


## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  data(trade, package = "fixest")
#  head(trade)

## ---- echo = FALSE, results = "asis"------------------------------------------
if(IS_FIXEST){
  data(trade, package = "fixest")
  tab = head(trade)
  knitr::kable(tab)
} else {
  msg()
}

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  plot_distr(trade$Euros)

## ---- echo = FALSE, eval = TRUE, fig.width=7----------------------------------
if(IS_FIXEST){
  plot_distr(trade$Euros)
} else {
  msg()
}

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  plot_distr(~ Euros | Origin, trade, mod.select = c("DE", "FR"))

## ---- echo = FALSE, eval = TRUE, fig.width=7----------------------------------
if(IS_FIXEST){
  plot_distr(~ Euros | Origin, trade, mod.select = c("DE", "FR"))
} else {
  msg()
}

## -----------------------------------------------------------------------------
plot_distr(iris$Petal.Length)

## -----------------------------------------------------------------------------
plot_distr(~ Petal.Length | Species, iris)

## -----------------------------------------------------------------------------
plot_distr(~ Petal.Length | Species, iris, mod.method = "split")

## -----------------------------------------------------------------------------
plot_distr(~ Sepal.Length | Species, iris, mod.method = "stack")

## -----------------------------------------------------------------------------
plot_distr(jnl_top_5p ~ institution, us_pub_econ, cumul = TRUE)

## ---- fig.width=7-------------------------------------------------------------
plot_lines(jnl_top_5p ~ year, us_pub_econ)

## ---- fig.width=7-------------------------------------------------------------
# Let's find out the top 3/5 institutions in terms of production
# we use plot_distr without plotting:
info = plot_distr(us_pub_econ$institution, plot = FALSE)
top3_instit = head(info$x, 3)
top5_instit = head(info$x, 5)
plot_lines(jnl_top_5p ~ year | institution, 
           us_pub_econ[institution %in% top3_instit])

## ---- fig.width=7-------------------------------------------------------------
plot_lines(1 ~ year | institution, us_pub_econ[institution %in% top5_instit])

## ---- fig.width=7-------------------------------------------------------------
plot_lines(1 ~ year | institution, us_pub_econ[institution %in% top5_instit], 
           smoothing_window = 1)

## -----------------------------------------------------------------------------
base_pub = us_pub_econ[, .(nb_pub_top_5 = sum(jnl_top_5p)), 
                          by = .(year, institution)]
plot_box(nb_pub_top_5 ~ year, base_pub)

## -----------------------------------------------------------------------------
base_pub[, isTop5 := institution %in% top5_instit]
plot_box(nb_pub_top_5 ~ year | isTop5, base_pub)

## -----------------------------------------------------------------------------
# we also drop the mean
plot_box(nb_pub_top_5 ~ year | isTop5, base_pub, outlier = FALSE)

## -----------------------------------------------------------------------------
plot_box(. ~ Species, iris)

