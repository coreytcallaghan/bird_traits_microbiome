# summarize the data

# packages
library(ggplot2)
library(GGally)
library(tidyverse)
library(forcats)
library(scales)

# prepare trait data
analysis_dat <- readRDS("Clean data/analysis_data.RDS") %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

length(unique(analysis_dat$ebird_COMMON_NAME))
nrow(analysis_dat)
length(unique(analysis_dat$DOI))

traits <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Mass, Range.Size, 
                habitat_breadth, Trophic.Niche, Trophic.Level, Habitat, 
                mean_flock_size, Migration) %>%
  distinct() %>%
  mutate(Mass_log10=log10(Mass)) %>%
  mutate(Range.Size_log10=log10(Range.Size)) %>%
  mutate(habitat_breadth_log10=log10(habitat_breadth)) %>%
  mutate(mean_flock_size_log10=log10(mean_flock_size))

length(unique(traits$ebird_COMMON_NAME))

traits %>%
  dplyr::select(Mass_log10, Range.Size_log10, habitat_breadth_log10, mean_flock_size_log10) %>%
  ggpairs(.)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

ggsave("Figures/continuous_predictors_ggpairs.png", width=7.1, height=7.1, units="in")

table(traits$Migration)
table(traits$Habitat)
table(traits$Trophic.Niche)
table(traits$Trophic.Level)

# make figure of number of species per family
analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, family) %>%
  distinct() %>%
  group_by(family) %>%
  summarize(N=n()) %>%
  arrange(N) %>%
  ggplot(., aes(x=fct_inorder(family), y=N))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Number of species")+
  xlab("")+
  scale_y_continuous(breaks=pretty_breaks())

ggsave("Figures/number_sp_per_family.png", width=7.3, height=7.1, units="in")
