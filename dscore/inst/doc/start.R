## ----setup, include=FALSE-----------------------------------------------------
  knitr::opts_chunk$set(echo = TRUE, fig.retina = 2)
  # options(width = 80)
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))

## ----graphkey, fig.width = 7, fig.height = 5, echo = FALSE--------------------
library(dscore)
ib <- builtin_itembank %>% 
  filter(key == "gsed") %>% 
  mutate(a = get_age_equivalent(items = item, pct = 50, itembank = builtin_itembank)$a * 12) %>% 
  select(a, instrument, label) %>% 
  na.omit()
  
ggplot(ib, aes(x = a, y = instrument, group = instrument)) +
  scale_y_discrete(limits = rev(unique(ib$instrument)), name = "") +
  scale_x_continuous(limits = c(0, 60), 
                     breaks = seq(0, 60, 12), name = "Age (months)") +
  geom_point(pch = 3, size = 1, colour = "blue") + 
  theme_light() +
  theme(axis.text.y = element_text(hjust = 0, family = "mono"))

## ----getlabels----------------------------------------------------------------
library(dscore)
get_labels("by3cgd018")

## ----decompose_itemnames------------------------------------------------------
decompose_itemnames(c("by3cgd018", "denfmd014"))

## ----table--------------------------------------------------------------------
items <- get_itemnames()
din <- decompose_itemnames(items)
knitr::kable(with(din, table(instrument, domain)), format = "html") %>% 
  kableExtra::column_spec(1, monospace = TRUE)

## ----ddigm--------------------------------------------------------------------
items <- head(get_itemnames(instrument = "mdt", domain = "gm"), 3)
get_labels(items)

## ----smalldataset-------------------------------------------------------------
data <- data.frame(id = c(1, 1, 2), age = c(1, 1.6, 0.9), mot1 = c(1, NA, NA), 
                   mot2 = c(0, 1, 1), mot3 = c(NA, 0, 1))
data

## ----rename-------------------------------------------------------------------
old_names <- names(data)[3:5]
new_names <- get_itemnames(instrument = "mdt", domain = "gm")[1:3]
names(data)[3:5] <- new_names
data

## ----milestones---------------------------------------------------------------
head(milestones[, c(1, 3, 4, 9:14)])

## -----------------------------------------------------------------------------
ds <- dscore(milestones)
dim(ds)

## -----------------------------------------------------------------------------
head(ds)

## ----bind---------------------------------------------------------------------
md <- cbind(milestones, ds)

## ----graphD, fig.width = 7, fig.height = 7------------------------------------
library(ggplot2)
library(dplyr)

r <- builtin_references %>% 
  filter(pop == "dutch") %>% 
  select(age, SDM2, SD0, SDP2)

ggplot(md, aes(x = a, y = d, group = id, color = sex)) + 
  theme_light() + 
  theme(legend.position = c(.85, .15)) +
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank()) +
  annotate("polygon", x = c(r$age, rev(r$age)), 
           y = c(r$SDM2, rev(r$SDP2)), alpha = 0.1, fill = "green") +
  annotate("line", x = r$age, y = r$SDM2, lwd = 0.3, alpha = 0.2, color = "green") +
  annotate("line", x = r$age, y = r$SDP2, lwd = 0.3, alpha = 0.2, color = "green") +
  annotate("line", x = r$age, y = r$SD0, lwd = 0.5, alpha = 0.2, color = "green") +
  coord_cartesian(xlim = c(0, 2.5)) +
  ylab(expression(paste(italic(D), "-score", sep = ""))) +
  xlab("Age (in years)") +
  scale_color_brewer(palette = "Set1") +
  geom_line(lwd = 0.1) +
  geom_point(size = 1)

## ----graphDAZ, fig.width = 7, fig.height = 5----------------------------------
ggplot(md, aes(x = a, y = daz, group = id, color = factor(sex))) + 
  theme_light() + 
  theme(legend.position = c(.85, .15)) +
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2, alpha = 0.1,
            fill = "green") +
  coord_cartesian(xlim = c(0, 2.5), 
                  ylim = c(-4, 4)) +
  geom_line(lwd = 0.1) +
  geom_point(size = 1) +
  xlab("Age (in years)") +
  ylab("DAZ") 

## ----density, fig.width = 7, fig.height = 4-----------------------------------
ggplot(md) + 
  theme_light() +
  scale_fill_brewer(palette = "Set1") +
  geom_density(aes(x = daz, group = sex, fill = sex), alpha = 0.4, 
               adjust = 1.5, color = "transparent") +
  xlab("DAZ")

## ----independence-------------------------------------------------------------
summary(lm(daz ~  age * sex, data = md))

## ----multilevel---------------------------------------------------------------
library(lme4)
lmer(daz ~  1 + age + sex + sex * age + (1 + age | id), data = md)

