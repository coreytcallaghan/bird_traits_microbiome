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

# a function that takes a predictor variable and fits a model
single_regression_model_function <- function(predictor_name){
  
  message(paste0("Modelling ", predictor_name))
  
  temp_dat_with_big_birds <- analysis_dat %>%
    dplyr::select(ebird_COMMON_NAME, gamma_Richness, DOI, predictor_name) %>%
    rename(predictor=4) %>%
    dplyr::filter(complete.cases(.))
  
  if(predictor_name %in% c("Mass", "Range.Size", "mean_flock_size", "habitat_breadth", "pop_abund")) {
    
    mod_with <- brms::brm(gamma_Richness ~ log10(predictor) + (1|DOI),
                          family=gaussian(),
                          data=temp_dat_with_big_birds,
                          warmup=1000,
                          iter=4000, 
                          chains=4,
                          control=list(adapt_delta=0.99),
                          file=paste0("model_objects/gamma_single_regression_richness_", predictor_name, ".RDS"))
    
    summary(mod_with)
    
    # coefficients table
    coefficients <- brms_SummaryTable(mod_with)
    
    # get a sample of draws from the posterior
    draws <- mod_with %>%
      spread_draws(b_log10predictor) %>%
      dplyr::select(b_log10predictor)
    
    saveRDS(draws, paste0("intermediate_results/gamma_single_regression_richness_draws_", predictor_name, ".RDS"))
    
    # compute R2
    r2 <- bayes_R2(mod_with, summary=TRUE) %>%
      as.data.frame()
    
    # get fixed effects plot
    fe_only <- tibble(predictor = seq(min(mod_with$data$predictor), 
                                      max(mod_with$data$predictor), length.out=100)) %>%
      add_fitted_draws(mod_with,
                       re_formula = NA,
                       scale = "response", n = 1e3)
    
    saveRDS(fe_only, paste0("intermediate_results/gamma_single_regression_richness_fe_only_", predictor_name, ".RDS"))
    
    fe_only_mean <- fe_only %>% 
      group_by(predictor) %>%
      summarize(.value = mean(.value))
    
    saveRDS(fe_only_mean, paste0("intermediate_results/gamma_single_regression_richness_fe_only_mean_", predictor_name, ".RDS"))
    
    out_df <- coefficients %>%
      mutate(model_R2=r2$Estimate) %>%
      mutate(number_species=length(unique(temp_dat_with_big_birds$ebird_COMMON_NAME))) %>%
      mutate(number_obs=nrow(temp_dat_with_big_birds)) %>%
      mutate(number_studies=length(unique(temp_dat_with_big_birds$DOI))) %>%
      mutate(type="with_big_birds") %>%
      mutate(predictor=predictor_name)
    
    saveRDS(out_df, paste0("intermediate_results/single_regression_richness_summary_", predictor_name, ".RDS"))
    
  } else if(predictor_name %in% c("Habitat")){
    
    mod_with <- brms::brm(gamma_Richness ~ predictor + (1|DOI),
                          family=gaussian(),
                          data=temp_dat_with_big_birds,
                          warmup=1000,
                          iter=4000, 
                          chains=4,
                          control=list(adapt_delta=0.99),
                          file=paste0("model_objects/gamma_single_regression_richness_", predictor_name, ".RDS"))
    
    summary(mod_with)
    
    # coefficients table
    coefficients <- brms_SummaryTable(mod_with)
    
    # draws <- mod_with %>%
    #   spread_draws(b_Intercept, b_predictorGrassland, b_predictorHumanModified, 
    #                b_predictorShrubland, b_predictorWetland, b_predictorWoodland)
    # 
    # saveRDS(draws, paste0("intermediate_results/gamma_single_regression_richness_draws_", predictor_name, ".RDS"))
    
    # compute R2
    r2 <- bayes_R2(mod_with, summary=TRUE) %>%
      as.data.frame()
    
    # get fixed effects plot
    fe_only <- tibble(predictor = rep_len(unique(mod_with$data$predictor), length.out=100)) %>%
      add_fitted_draws(mod_with,
                       re_formula = NA,
                       scale = "response", n = 1e3)
    
    saveRDS(fe_only, paste0("intermediate_results/gamma_single_regression_richness_fe_only_", predictor_name, ".RDS"))
    
    fe_only_mean <- fe_only %>% 
      group_by(predictor) %>%
      summarize(.value = mean(.value))
    
    saveRDS(fe_only_mean, paste0("intermediate_results/gamma_single_regression_richness_fe_only_mean_", predictor_name, ".RDS"))
    
    out_df <- coefficients %>%
      mutate(model_R2=r2$Estimate) %>%
      mutate(number_species=length(unique(temp_dat_with_big_birds$ebird_COMMON_NAME))) %>%
      mutate(number_obs=nrow(temp_dat_with_big_birds)) %>%
      mutate(number_studies=length(unique(temp_dat_with_big_birds$DOI))) %>%
      mutate(type="with_big_birds") %>%
      mutate(predictor=predictor_name)
    
    saveRDS(out_df, paste0("intermediate_results/gamma_single_regression_richness_summary_", predictor_name, ".RDS"))
    
  } else if(predictor_name %in% c("Primary.Lifestyle")){
    
    mod_with <- brms::brm(gamma_Richness ~ predictor + (1|DOI),
                          family=gaussian(),
                          data=temp_dat_with_big_birds,
                          warmup=1000,
                          iter=4000, 
                          chains=4,
                          control=list(adapt_delta=0.99),
                          file=paste0("model_objects/gamma_single_regression_richness_", predictor_name, ".RDS"))
    
    summary(mod_with)
    
    # coefficients table
    coefficients <- brms_SummaryTable(mod_with)
    
    # draws <- mod_with %>%
    #   spread_draws(b_Intercept, b_predictorTerrestrial, b_predictorInsessorial,
    #                b_predictorGeneralist)
    # 
    # saveRDS(draws, paste0("intermediate_results/gamma_single_regression_richness_draws_", predictor_name, ".RDS"))
    
    # compute R2
    r2 <- bayes_R2(mod_with, summary=TRUE) %>%
      as.data.frame()
    
    # get fixed effects plot
    fe_only <- tibble(predictor = rep_len(unique(mod_with$data$predictor), length.out=100)) %>%
      add_fitted_draws(mod_with,
                       re_formula = NA,
                       scale = "response", n = 1e3)
    
    saveRDS(fe_only, paste0("intermediate_results/gamma_single_regression_richness_fe_only_", predictor_name, ".RDS"))
    
    fe_only_mean <- fe_only %>% 
      group_by(predictor) %>%
      summarize(.value = mean(.value))
    
    saveRDS(fe_only_mean, paste0("intermediate_results/gamma_single_regression_richness_fe_only_mean_", predictor_name, ".RDS"))
    
    out_df <- coefficients %>%
      mutate(model_R2=r2$Estimate) %>%
      mutate(number_species=length(unique(temp_dat_with_big_birds$ebird_COMMON_NAME))) %>%
      mutate(number_obs=nrow(temp_dat_with_big_birds)) %>%
      mutate(number_studies=length(unique(temp_dat_with_big_birds$DOI))) %>%
      mutate(type="with_big_birds") %>%
      mutate(predictor=predictor_name)
    
    saveRDS(out_df, paste0("intermediate_results/gamma_single_regression_richness_summary_", predictor_name, ".RDS"))
    
  } else if(predictor_name %in% c("Trophic.Niche")){
    
    mod_with <- brms::brm(gamma_Richness ~ predictor + (1|DOI),
                          family=gaussian(),
                          data=temp_dat_with_big_birds,
                          warmup=1000,
                          iter=4000, 
                          chains=4,
                          control=list(adapt_delta=0.99),
                          file=paste0("model_objects/gamma_single_regression_richness_", predictor_name, ".RDS"))
    
    summary(mod_with)
    
    # coefficients table
    coefficients <- brms_SummaryTable(mod_with)
    
    # draws <- mod_with %>%
    #   spread_draws(b_Intercept, b_predictorGranivore, b_predictorOmnivore,
    #                b_predictorInvertivore, b_predictorFrugivore, b_predictorNectarivore, 
    #                b_predictorHerbivoreterrestrial)
    # 
    # saveRDS(draws, paste0("intermediate_results/gamma_single_regression_richness_draws_", predictor_name, ".RDS"))
    # 
    # compute R2
    r2 <- bayes_R2(mod_with, summary=TRUE) %>%
      as.data.frame()
    
    # get fixed effects plot
    fe_only <- tibble(predictor = rep_len(unique(mod_with$data$predictor), length.out=100)) %>%
      add_fitted_draws(mod_with,
                       re_formula = NA,
                       scale = "response", n = 1e3)
    
    saveRDS(fe_only, paste0("intermediate_results/gamma_single_regression_richness_fe_only_", predictor_name, ".RDS"))
    
    fe_only_mean <- fe_only %>% 
      group_by(predictor) %>%
      summarize(.value = mean(.value))
    
    saveRDS(fe_only_mean, paste0("intermediate_results/gamma_single_regression_richness_fe_only_mean_", predictor_name, ".RDS"))
    
    out_df <- coefficients %>%
      mutate(model_R2=r2$Estimate) %>%
      mutate(number_species=length(unique(temp_dat_with_big_birds$ebird_COMMON_NAME))) %>%
      mutate(number_obs=nrow(temp_dat_with_big_birds)) %>%
      mutate(number_studies=length(unique(temp_dat_with_big_birds$DOI))) %>%
      mutate(type="with_big_birds") %>%
      mutate(predictor=predictor_name)
    
    saveRDS(out_df, paste0("intermediate_results/gamma_single_regression_richness_summary_", predictor_name, ".RDS"))
    
  } else if(predictor_name %in% c("Trophic.Level")){
    
    mod_with <- brms::brm(gamma_Richness ~ predictor + (1|DOI),
                          family=gaussian(),
                          data=temp_dat_with_big_birds,
                          warmup=1000,
                          iter=4000, 
                          chains=4,
                          control=list(adapt_delta=0.99),
                          file=paste0("model_objects/gamma_single_regression_richness_", predictor_name, ".RDS"))
    
    summary(mod_with)
    
    # coefficients table
    coefficients <- brms_SummaryTable(mod_with)
    
    draws <- mod_with %>%
      spread_draws(b_Intercept, b_predictorOmnivore, b_predictorHerbivore)
    
    saveRDS(draws, paste0("intermediate_results/gamma_single_regression_richness_draws_", predictor_name, ".RDS"))
    
    # compute R2
    r2 <- bayes_R2(mod_with, summary=TRUE) %>%
      as.data.frame()
    
    # get fixed effects plot
    fe_only <- tibble(predictor = rep_len(unique(mod_with$data$predictor), length.out=100)) %>%
      add_fitted_draws(mod_with,
                       re_formula = NA,
                       scale = "response", n = 1e3)
    
    saveRDS(fe_only, paste0("intermediate_results/gamma_single_regression_richness_fe_only_", predictor_name, ".RDS"))
    
    fe_only_mean <- fe_only %>% 
      group_by(predictor) %>%
      summarize(.value = mean(.value))
    
    saveRDS(fe_only_mean, paste0("intermediate_results/gamma_single_regression_richness_fe_only_mean_", predictor_name, ".RDS"))
    
    out_df <- coefficients %>%
      mutate(model_R2=r2$Estimate) %>%
      mutate(number_species=length(unique(temp_dat_with_big_birds$ebird_COMMON_NAME))) %>%
      mutate(number_obs=nrow(temp_dat_with_big_birds)) %>%
      mutate(number_studies=length(unique(temp_dat_with_big_birds$DOI))) %>%
      mutate(type="with_big_birds") %>%
      mutate(predictor=predictor_name)
    
    saveRDS(out_df, paste0("intermediate_results/gamma_single_regression_richness_summary_", predictor_name, ".RDS"))
    
  }
  
}

# now fit a model for every trait we are interested in
lapply(c("Mass", "Range.Size", "habitat_breadth", "pop_abund", 
         "mean_flock_size", "Habitat", "Trophic.Level", "Trophic.Niche",
         "Primary.Lifestyle"), single_regression_model_function)








