# make a table of the sample sizes for each individual model
# packages
library(dplyr)
library(ggplot2)
library(ggstance)
library(tidybayes)
library(tidyr)
library(forcats)
library(patchwork)

# prepare trait data
analysis_dat <- readRDS("Clean data/analysis_data_alpha.RDS") %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

length(unique(analysis_dat$ebird_COMMON_NAME))
nrow(analysis_dat)
length(unique(analysis_dat$DOI))


# first get a summary of all traits
get_summary <- function(predictor_name){
  
  dat <- readRDS(paste0("intermediate_results/single_regression_richness_summary_", predictor_name, ".RDS"))
  
  return(dat)
  
}

overall_summary <- bind_rows(lapply(c("Mass", "Range.Size", "habitat_breadth", "pop_abund", 
                                      "mean_flock_size", "Habitat", "Trophic.Level", 
                                      "Primary.Lifestyle"), get_summary))

table <- overall_summary %>%
  mutate(predictor2=case_when(predictor=="Mass" ~ "Body mass",
                              predictor=="mean_flock_size" ~ "Flock size",
                              predictor=="Range.Size" ~ "Range size",
                              predictor=="Habitat" ~ "Primary habitat",
                              predictor=="Trophic.Level" ~ "Trophic level",
                              predictor=="Primary.Lifestyle" ~ "Primary lifestyle",
                              predictor=="pop_abund" ~ "Global abundance",
                              predictor=="habitat_breadth" ~ "Habitat breadth")) %>%
  dplyr::select(predictor2, number_species, number_obs, number_studies) %>%
  distinct()






######################################
######################################
########### GAMMA #####################
######################################
######################################
# prepare trait data
# prepare trait data
analysis_dat <- readRDS("Clean data/analysis_data_gamma.RDS") %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

length(unique(analysis_dat$ebird_COMMON_NAME))
nrow(analysis_dat)
length(unique(analysis_dat$DOI))


# first get a summary of all traits
get_summary <- function(predictor_name){
  
  dat <- readRDS(paste0("intermediate_results/gamma_single_regression_richness_summary_", predictor_name, ".RDS"))
  
  return(dat)
  
}

overall_summary <- bind_rows(lapply(c("Mass", "Range.Size", "habitat_breadth", "pop_abund", 
                                      "mean_flock_size", "Habitat", "Trophic.Level",
                                      "Primary.Lifestyle"), get_summary))

table_gamma <- overall_summary %>%
  mutate(predictor2=case_when(predictor=="Mass" ~ "Body mass",
                              predictor=="mean_flock_size" ~ "Flock size",
                              predictor=="Range.Size" ~ "Range size",
                              predictor=="Habitat" ~ "Primary habitat",
                              predictor=="Trophic.Level" ~ "Trophic level",
                              predictor=="Primary.Lifestyle" ~ "Primary lifestyle",
                              predictor=="pop_abund" ~ "Global abundance",
                              predictor=="habitat_breadth" ~ "Habitat breadth")) %>%
  dplyr::select(predictor2, number_species, number_obs, number_studies) %>%
  distinct()

