# summarize multiple regression outputs

library(tidyverse)
library(forcats)
library(ggrepel)


# first get a summary of all traits
get_summary <- function(predictor_name){
  
  dat <- readRDS(paste0("intermediate_results/single_regression_ISimpson_summary_", predictor_name, ".RDS"))
  
  return(dat)
  
}

overall_summary <- bind_rows(lapply(c("Mass", "Range.Size", "habitat_breadth", "pop_abund", 
                                      "mean_flock_size", "Habitat", "Trophic.Level",
                                      "Primary.Lifestyle"), get_summary))

summary_simpson <- overall_summary %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
  mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
  mutate(predictor2=case_when(predictor=="Mass" ~ "Body mass",
                              predictor=="mean_flock_size" ~ "Flock size",
                              predictor=="Range.Size" ~ "Range size",
                              predictor=="Habitat" ~ "Primary habitat",
                              predictor=="Trophic.Level" ~ "Trophic level",
                              predictor=="Primary.Lifestyle" ~ "Primary lifestyle",
                              predictor=="pop_abund" ~ "Global abundance",
                              predictor=="habitat_breadth" ~ "Habitat breadth")) %>%
  mutate(Covariate=gsub("predictor", "", Covariate)) %>%
  mutate(Covariate=gsub("HumanModified", "Human modified", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="Mass" & Covariate=="log10", "Body mass", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="mean_flock_size" & Covariate=="log10", "Flock size", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="pop_abund" & Covariate=="log10", "Global abundance", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="habitat_breadth" & Covariate=="log10", "Habitat breadth", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="Range.Size" & Covariate=="log10", "Range size", Covariate)) %>%
  mutate(predictor2=factor(predictor2, levels=c("Body mass", "Flock size", "Global abundance",
                                                "Habitat breadth", "Range size", "Trophic level",
                                                "Primary habitat", "Primary lifestyle"))) %>%
  arrange(Estimate) %>%
  mutate(analysis="Simpson")

# first get a summary of all traits
get_summary <- function(predictor_name){
  
  dat <- readRDS(paste0("intermediate_results/single_regression_richness_summary_", predictor_name, ".RDS"))
  
  return(dat)
  
}

overall_summary <- bind_rows(lapply(c("Mass", "Range.Size", "habitat_breadth", "pop_abund", 
                                      "mean_flock_size", "Habitat", "Trophic.Level", 
                                      "Primary.Lifestyle"), get_summary))

summary_richness <- overall_summary %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
  mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
  mutate(predictor2=case_when(predictor=="Mass" ~ "Body mass",
                              predictor=="mean_flock_size" ~ "Flock size",
                              predictor=="Range.Size" ~ "Range size",
                              predictor=="Habitat" ~ "Primary habitat",
                              predictor=="Trophic.Level" ~ "Trophic level",
                              predictor=="Primary.Lifestyle" ~ "Primary lifestyle",
                              predictor=="pop_abund" ~ "Global abundance",
                              predictor=="habitat_breadth" ~ "Habitat breadth")) %>%
  mutate(Covariate=gsub("predictor", "", Covariate)) %>%
  mutate(Covariate=gsub("HumanModified", "Human modified", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="Mass" & Covariate=="log10", "Body mass", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="mean_flock_size" & Covariate=="log10", "Flock size", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="pop_abund" & Covariate=="log10", "Global abundance", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="habitat_breadth" & Covariate=="log10", "Habitat breadth", Covariate)) %>%
  mutate(Covariate=ifelse(predictor=="Range.Size" & Covariate=="log10", "Range size", Covariate)) %>%
  mutate(predictor2=factor(predictor2, levels=c("Body mass", "Flock size", "Global abundance",
                                                "Habitat breadth", "Range size", "Trophic level",
                                                "Primary habitat", "Primary lifestyle"))) %>%
  arrange(Estimate) %>%
  mutate(analysis="Richness")



# plot the simpson vs richness values
summary_richness %>%
  dplyr::select(Covariate, predictor2, Estimate, upr, lwr) %>%
  rename(richness=Estimate,
         richness_upr=upr,
         richness_lwr=lwr) %>%
  left_join(., summary_simpson %>%
              dplyr::select(Covariate, predictor2, Estimate, upr, lwr) %>%
              rename(simpson=Estimate,
                     simpson_upr=upr,
                     simpson_lwr=lwr)) %>%
  dplyr::filter(Covariate != "Intercept") %>%
  ggplot(., aes(x=richness, y=simpson, label=Covariate))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Richness parameter estimate")+
  ylab("Simpson parameter estimate")+
  geom_hline(yintercept=0, color="red", linetype="dashed")+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  geom_smooth(method="lm")

ggsave("Figures/richness_vs_simpson_estimates_single_model_results.png", height=5.6, width=6.1, units="in")
