## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ChangePointTaylor)
library(dplyr)
library(ggplot2)

## ----paged.print=TRUE---------------------------------------------------------
US_Trade_Deficit

## ----fig.height=4, fig.width=7------------------------------------------------

trade_deficit_plot <- US_Trade_Deficit %>%
  mutate(date = as.Date(paste(date, "1"), format = "%b '%y %d")) %>%
  ggplot(aes(x = date, y = deficit_billions, group = 1)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b '%y") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust =1),
    axis.title.x = element_blank()
  ) +
  ggtitle("US Trade Deficit: 1987-1988")
  
trade_deficit_plot

## -----------------------------------------------------------------------------
change_point_analyzer(US_Trade_Deficit$deficit_billions)

## -----------------------------------------------------------------------------
change_points <- change_point_analyzer(US_Trade_Deficit$deficit_billions, label = US_Trade_Deficit$date)
change_points

## ----fig.height=4, fig.width=7------------------------------------------------
trade_deficit_plot +
  geom_vline(xintercept = as.Date(paste(change_points$label, "1"), format = "%b '%y %d"), color = "steelblue", linetype = "dashed", size = 1.3)

## ----message=FALSE, warning=FALSE---------------------------------------------
bench::mark(
  change_point_analyzer(US_Trade_Deficit$deficit_billions, label = US_Trade_Deficit$date, n_bootstraps = 1000)
 ,change_point_analyzer(US_Trade_Deficit$deficit_billions, label = US_Trade_Deficit$date, n_bootstraps = 10000)
  ,check = F
 ,min_iterations = 2
 ,max_iterations = 5
) %>%
  mutate(expression = c("1000 Bootstraps", "10000 Bootstraps")) %>%
  select(expression:mem_alloc)

## -----------------------------------------------------------------------------
change_point_analyzer(US_Trade_Deficit$deficit_billions, label = US_Trade_Deficit$date, min_candidate_conf = 0.66,  min_tbl_conf = 0.95)

