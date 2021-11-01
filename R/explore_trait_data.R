# explore trait data

library(dplyr)
library(ggplot2)
library(naniar)

traits <- readRDS("bird_trait_predictors.RDS")

vis_miss(traits)+
  coord_flip()
