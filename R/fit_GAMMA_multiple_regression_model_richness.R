# this script is to fit models with
# in brms
# as single regression models

source("R/global_functions.R")

# packages
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(lme4)
library(ggeffects)
library(naniar)
library(performance)
library(brms)
library(tidybayes)

# prepare trait data
analysis_dat <- readRDS("Clean data/analysis_data_gamma.RDS") %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

length(unique(analysis_dat$ebird_COMMON_NAME))
nrow(analysis_dat)
length(unique(analysis_dat$DOI))

model_dat <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, gamma_Richness, gamma_ISimpson, DOI, Sample_type,
                Mass, Range.Size, mean_flock_size, habitat_breadth,
                pop_abund, Habitat, Trophic.Level, Trophic.Niche,
                Primary.Lifestyle) %>%
  dplyr::filter(complete.cases(.))

length(unique(model_dat$ebird_COMMON_NAME))
nrow(model_dat)
length(unique(model_dat$DOI))

mod_with <- brms::brm(gamma_Richness ~ log10(Mass) + log10(Range.Size) + log10(habitat_breadth) +
                        log10(pop_abund) + log10(mean_flock_size) + Habitat + Trophic.Level +
                        Primary.Lifestyle + (1|DOI/Sample_type) + (1|ebird_COMMON_NAME),
                      family=gaussian(),
                      data=model_dat,
                      warmup=1000,
                      iter=4000, 
                      chains=4,
                      control=list(adapt_delta=0.99),
                      file=paste0("model_objects/gamma_multiple_regression_richness.RDS"))

summary(mod_with)

# coefficients table
coefficients <- brms_SummaryTable(mod_with)

# get a sample of draws from the posterior
draws <- mod_with %>%
  spread_draws(b_Intercept, 
               b_log10Mass, 
               b_log10Range.Size,
               b_log10habitat_breadth, 
               b_log10pop_abund,
               b_log10mean_flock_size, 
               b_HabitatGrassland, b_HabitatHumanModified, b_HabitatShrubland, b_HabitatWoodland, 
               b_Trophic.LevelHerbivore, b_Trophic.LevelOmnivore, 
               b_Primary.LifestyleInsessorial, b_Primary.LifestyleTerrestrial)

saveRDS(draws, paste0("intermediate_results/gamma_multiple_regression_richness_draws.RDS"))

# compute R2
r2 <- bayes_R2(mod_with, summary=TRUE) %>%
  as.data.frame()

out_df <- coefficients %>%
  mutate(model_R2=r2$Estimate) %>%
  mutate(number_species=length(unique(model_dat$ebird_COMMON_NAME))) %>%
  mutate(number_obs=nrow(model_dat)) %>%
  mutate(number_studies=length(unique(model_dat$DOI))) %>%
  mutate(type="with_big_birds")

saveRDS(out_df, paste0("intermediate_results/gamma_multiple_regression_richness_summary.RDS"))

# # get fixed effects plot
# fe_only <- tibble(habitat_breadth = seq(min(mod_with$data$habitat_breadth), 
#                                                  max(mod_with$data$habitat_breadth), length.out=100)) %>%
#   mutate(Mass=mean(mod_with$data$Mass)) %>%
#   mutate(Range.Size=mean(mod_with$data$Range.Size)) %>%
#   mutate(pop_abund=mean(mod_with$data$pop_abund)) %>%
#   mutate(mean_flock_size=mean(mod_with$data$mean_flock_size)) %>%
#   crossing(Habitat=unique(mod_with$data$Habitat)) %>%
#   crossing(Trophic.Niche=unique(mod_with$data$Trophic.Niche)) %>%
#   crossing(Trophic.Level=unique(mod_with$data$Trophic.Level)) %>%
#   crossing(Primary.Lifestyle=unique(mod_with$data$Primary.Lifestyle)) %>%
#   add_fitted_draws(mod_with,
#                    re_formula = NA,
#                    scale = "response", n = 100)
# 
# 
# saveRDS(fe_only, paste0("intermediate_results/single_regression_richness_fe_only_", predictor_name, ".RDS"))
# 
# fe_only_mean <- fe_only %>%
#   ungroup() %>%
#   group_by(habitat_breadth) %>%
#   summarize(sd=sd(.value),
#             mean = mean(.value),
#             N=n())
# 
# saveRDS(fe_only_mean, paste0("intermediate_results/single_regression_richness_fe_only_mean_", predictor_name, ".RDS"))
# 
# 
