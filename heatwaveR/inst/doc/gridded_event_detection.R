## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(fig.width = 4, fig.align = 'center',
                      echo = TRUE, warning = FALSE, message = FALSE, 
                      eval = FALSE, tidy = FALSE)

## ----load-pkg-----------------------------------------------------------------
#  library(dplyr) # For basic data manipulation
#  library(ggplot2) # For visualising data
#  library(heatwaveR) # For detecting MHWs
#  library(tidync) # For easily dealing with NetCDF data
#  library(doParallel) # For parallel processing

## ----load-data----------------------------------------------------------------
#  OISST <- readRDS("~/Desktop/OISST_vignette.Rds")

## ----detect-func--------------------------------------------------------------
#  event_only <- function(df){
#    # First calculate the climatologies
#    clim <- ts2clm(data = df, climatologyPeriod = c("1982-01-01", "2011-01-01"))
#    # Then the events
#    event <- detect_event(data = clim)
#    # Return only the event metric dataframe of results
#    return(event$event)
#  }

## ----detect-purrr-------------------------------------------------------------
#  system.time(
#  # First we start by choosing the 'OISST' dataframe
#  MHW_dplyr <- OISST %>%
#    # Then we group the data by the 'lon' and 'lat' columns
#    group_by(lon, lat) %>%
#    # Then we run our MHW detecting function on each group
#    group_modify(~event_only(.x))
#  ) # ~123 seconds

## ----detect-plyr--------------------------------------------------------------
#  # NB: One should never use ALL available cores, save at least 1 for other essential tasks
#  # The computer I'm writing this vignette on has 8 cores, so I use 7 here
#  registerDoParallel(cores = 7)
#  
#  # Detect events
#  system.time(
#  MHW_plyr <- plyr::ddply(.data = OISST, .variables = c("lon", "lat"), .fun = event_only, .parallel = TRUE)
#  ) # 33 seconds

## ----lon-files----------------------------------------------------------------
#  for(i in 1:length(unique(OISST$lon))){
#    OISST_sub <- OISST %>%
#      filter(lon == unique(lon)[i])
#    saveRDS(object = OISST_sub, file = paste0("~/Desktop/OISST_lon_",i,".Rds"))
#  }

## ----detect-both--------------------------------------------------------------
#  # The 'dplyr' wrapper function to pass to 'plyr'
#  dplyr_wraper <- function(file_name){
#    MHW_dplyr <- readRDS(file_name) %>%
#      group_by(lon, lat) %>%
#      group_modify(~event_only(.x))
#  }
#  # Create a vector of the files we want to use
#  OISST_files <- dir("~/Desktop", pattern = "OISST_lon_*", full.names = T)
#  
#  # Use 'plyr' technique to run 'dplyr' technique with multiple cores
#  system.time(
#  MHW_result <- plyr::ldply(OISST_files, .fun = dplyr_wraper, .parallel = T)
#  ) # 31 seconds

## ----event-tally--------------------------------------------------------------
#  # summarise the number of unique longitude, latitude and year combination:
#  event_freq <- MHW_result %>%
#    mutate(year = lubridate::year(date_start)) %>%
#    group_by(lon, lat, year) %>%
#    summarise(n = n(), .groups = "drop")
#  head(event_freq)
#  
#  # create complete grid for merging with:
#  sst_grid <- OISST %>%
#    select(lon, lat, t) %>%
#    mutate(t = lubridate::year(t)) %>%
#    dplyr::rename(year = t) %>%
#    distinct()
#  
#  # and merge:
#  OISST_n <- left_join(sst_grid, event_freq, by = c("lon", "lat", "year")) %>%
#    mutate(n = ifelse(is.na(n), 0, n))

## ----trend-fun----------------------------------------------------------------
#  lin_fun <- function(ev) {
#    mod1 <- glm(n ~ year, family = poisson(link = "log"), data = ev)
#    # extract slope coefficient and its p-value
#    tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
#                     p = summary(mod1)$coefficients[2,4])
#    return(tr)
#  }

## ----apply-trend-fun-plyr-----------------------------------------------------
#  OISST_nTrend <- plyr::ddply(OISST_n, c("lon", "lat"), lin_fun, .parallel = T)
#  OISST_nTrend$pval <- cut(OISST_nTrend$p, breaks = c(0, 0.001, 0.01, 0.05, 1))
#  head(OISST_nTrend)

## ----prep-geography-----------------------------------------------------------
#  # The base map
#  map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
#    dplyr::rename(lon = long)

## ----the-figure---------------------------------------------------------------
#  map_slope <- ggplot(OISST_nTrend, aes(x = lon, y = lat)) +
#    geom_rect(size = 0.2, fill = NA,
#         aes(xmin = lon - 0.1, xmax = lon + 0.1, ymin = lat - 0.1, ymax = lat + 0.1,
#             colour = pval)) +
#    geom_raster(aes(fill = slope), interpolate = FALSE, alpha = 0.9) +
#    scale_fill_gradient2(name = "count/year (slope)", high = "red", mid = "white",
#                         low = "darkblue", midpoint = 0,
#                         guide = guide_colourbar(direction = "horizontal",
#                                                 title.position = "top")) +
#    scale_colour_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,1]"),
#                        values = c("firebrick1", "firebrick2", "firebrick3", "white"),
#                        name = "p-value", guide = FALSE) +
#    geom_polygon(data = map_base, aes(group = group),
#                 colour = NA, fill = "grey80") +
#    coord_fixed(ratio = 1, xlim = c(13.0, 23.0), ylim = c(-33, -42), expand = TRUE) +
#    labs(x = "", y = "") +
#    theme_bw() +
#    theme(legend.position = "bottom")
#  
#  map_p <- ggplot(OISST_nTrend, aes(x = lon, y = lat)) +
#    geom_raster(aes(fill = pval), interpolate = FALSE) +
#    scale_fill_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]",
#                                 "(0.05,0.1]", "(0.1,0.5]", "(0.5,1]"),
#                      values = c("black", "grey20", "grey40",
#                                 "grey80", "grey90", "white"),
#                      name = "p-value",
#                      guide = guide_legend(direction = "horizontal",
#                                                 title.position = "top")) +
#    geom_polygon(data = map_base, aes(group = group),
#                 colour = NA, fill = "grey80") +
#    coord_fixed(ratio = 1, xlim = c(13.0, 23.0), ylim = c(-33, -42), expand = TRUE) +
#    labs(x = "", y = "") +
#    theme_bw() +
#    theme(legend.position = "bottom")
#  
#  map_both <- ggpubr::ggarrange(map_slope, map_p, align = "hv")
#  map_both

