## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(fig.width = 8, fig.height = 3, fig.align = 'centre',
                      echo = TRUE, warning = FALSE, message = FALSE,
                      eval = TRUE, tidy = FALSE)

## ----load-libs----------------------------------------------------------------
library(dplyr)
library(ggpubr)
library(heatwaveR)

## ----data-prep, eval=T--------------------------------------------------------
Algiers <- Algiers

## ----clim-calc----------------------------------------------------------------
# The tMax threshold
# The current WMO standard climatology period is 1981-01-01 to 2010-12-31 and should be used when possible
# We rather use 1961-01-01 to 1990-01-01 as this is the oldest 30 year period available in the data
tMax_clim <- ts2clm(data = Algiers, y = tMax, climatologyPeriod = c("1961-01-01", "1990-12-31"), pctile = 90)

# The tMin exceedance
# Note the use here of 'minDuration = 3' and 'maxGap = 1' as the default atmospheric arguments
# The default marine arguments are 'minDuration = 5' and 'maxGap = 2'
tMin_exc <- exceedance(data = Algiers, y = tMin, threshold = 19, minDuration = 3, maxGap = 1)$threshold

## ----events-------------------------------------------------------------------
# Note that because we calculated our 90th percentile threshold on a column named 'tMax' 
# and not the default column name 'temp', we must specify this below with 'y = tMax'
events <- detect_event(data = tMax_clim, y = tMax, # The 90th percentile threshold
                       threshClim2 = tMin_exc$exceedance) # The flat exceedance threshold

## ----visuals, fig.height=6----------------------------------------------------
# The code to create a bubble plot for the heatwave results
bubble_plot <- ggplot(data = events$event, aes(x = date_peak, y = intensity_max)) +
  geom_point(aes(size = intensity_cumulative), shape = 21, fill = "salmon", alpha = 0.8) +
  labs(x = NULL, y = "Maximum Intensity [°C] ", size = "Cumulative Intensity [°C x days]") +
  scale_size_continuous(range = c(1, 10), 
                        guide = guide_legend(title.position = "top", direction = "horizontal")) +
  theme_bw() +
  theme(legend.position = c(0.3, 0.12),
        legend.box.background = element_rect(colour = "black"))

# Don't forget to set 'event_line(y = tMax)'
ggarrange(event_line(events, y = tMax, metric = "intensity_max"),
          event_line(events, y = tMax, metric = "intensity_max", category = T),
          lolli_plot(events),
          bubble_plot,
          ncol = 2, nrow = 2, align = "hv")

## ----alt-two-thresh-calc------------------------------------------------------
# Note that because we are not using the standard column name 'temp' we must
# specify the chosen column name twice, once for ts2clm() and again for detect_event()

# First threshold based on tMin
thresh_tMin <- ts2clm(data = Algiers, y = tMin, pctile = 80, 
                      climatologyPeriod = c("1961-01-01", "1990-12-31"))

# Second threshold based on tMax
# Be careful here that you put the arguments within the correct brackets
thresh_tMax <- detect_event(ts2clm(data = Algiers, y = tMax, pctile = 90, 
                                   climatologyPeriod = c("1961-01-01", "1990-12-31")),
                            # These arguments are passed to detect_event(), not ts2clm()
                            minDuration = 3, maxGap = 0, y = tMax, protoEvents = T)

# Detect/calculate events using the two precalculated thresholds
# Because detect_event() is not able to deduce which arguments we used above,
# we must again tell it explicitly here
events_two_thresh <- detect_event(data = thresh_tMin, y = tMin, minDuration = 10, maxGap = 2,
                                  threshClim2 = thresh_tMax$event, minDuration2 = 3, maxGap2 = 0)

# Or to simply use one threshold
events_one_thresh <- detect_event(data = thresh_tMin, y = tMin, minDuration = 10, maxGap = 2)

## ----alt-two-thresh-lollis----------------------------------------------------
ggarrange(lolli_plot(events_one_thresh), lolli_plot(events_two_thresh), labels = c("One threshold", "Two thresholds"))

## ----alt-two-thresh-events----------------------------------------------------
head(events_one_thresh$event)
head(events_two_thresh$event)

## ----event-data-frame---------------------------------------------------------
# Pull out each data.frame as their own object for easier use
events_one_event <- events_one_thresh$event
events_one_climatology <- events_one_thresh$climatology

## ----filter-join--------------------------------------------------------------
# Join the two threshold dataframes
two_thresh <- left_join(events_one_climatology, thresh_tMax, by = c("t"))

# Remove all days that did not qualify as events in both thresholds
two_thresh_filtered <- two_thresh %>%
  filter(event.x == TRUE,
         event.y == TRUE)

## ----filter-one-thresh--------------------------------------------------------
# Copy data with a new name
events_one_thresh_filtered <- events_one_thresh

# Then filter
events_one_thresh_filtered$event <- events_one_thresh_filtered$event %>% 
  filter(event_no %in% two_thresh_filtered$event_no.x)

# Compare results
head(events_one_thresh_filtered$event)
head(events_two_thresh$event)

## ----lolliplot-duration, fig.cap="Difference in duration (days) of events given different applications of thresholds. Note the difference in the y-axes."----
ggarrange(lolli_plot(events_two_thresh, metric = "duration"), 
          lolli_plot(events_one_thresh_filtered, metric = "duration"), 
          labels = c("Double threshold", "Filter threshold"))

## ----lolliplot-int-cum, fig.cap="Difference in cumulative intensity (°C x days) of events given different applications of thresholds. Note the difference in the y-axes."----
ggarrange(lolli_plot(events_two_thresh, metric = "intensity_cumulative"), 
          lolli_plot(events_one_thresh_filtered, metric = "intensity_cumulative"), 
          labels = c("Double threshold", "Filter threshold"))

