# summarize multiple regression outputs

library(tidyverse)
library(forcats)
library(ggrepel)

mult_results_richness <- readRDS("intermediate_results/multiple_regression_richness_summary.RDS") %>%
  mutate(scale="alpha") %>%
  mutate(response="richness")

mult_results_simpson <- readRDS("intermediate_results/multiple_regression_simpson_summary.RDS") %>%
  mutate(scale="alpha") %>%
  mutate(response="simpson")


mult_results_richness %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
  mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
  mutate(Covariate=gsub("log10", "", Covariate)) %>%
  mutate(Covariate=gsub("Habitat", "Habitat - ", Covariate)) %>%
  mutate(Covariate=gsub("mean_flock_size", "Flock size", Covariate)) %>%
  mutate(Covariate=gsub("pop_abund", "Global abundance", Covariate)) %>%
  mutate(Covariate=gsub("HumanModified", "Human modified", Covariate)) %>%
  mutate(Covariate=gsub("Primary.Lifestyle", "Primary lifestyle - ", Covariate)) %>%
  mutate(Covariate=gsub("Range.Size", "Range size", Covariate)) %>%
  mutate(Covariate=gsub("habitat_breadth", "Habitat breadth", Covariate)) %>%
  mutate(Covariate=gsub("Trophic.Level", "Trophic level - ", Covariate)) %>%
  mutate(Covariate=gsub("Mass", "Body mass", Covariate)) %>%
  arrange(Estimate) %>%
  ggplot(., aes(x=fct_inorder(Covariate), y=Estimate))+
  geom_point()+
  geom_errorbar(aes(x=fct_inorder(Covariate), y=Estimate, ymin=lwr, ymax=upr), width=0)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red", linetype="dashed")+
  coord_flip()+
  ylab("Parameter estimate")+
  xlab("")+
  ggtitle("Alpha Richness")

ggsave("Figures/alpha_richness_multiple_regression_results.png", width=6.9, height=6.4, units="in")
  
mult_results_simpson %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
  mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
  mutate(Covariate=gsub("log10", "", Covariate)) %>%
  mutate(Covariate=gsub("Habitat", "Habitat - ", Covariate)) %>%
  mutate(Covariate=gsub("mean_flock_size", "Flock size", Covariate)) %>%
  mutate(Covariate=gsub("pop_abund", "Global abundance", Covariate)) %>%
  mutate(Covariate=gsub("HumanModified", "Human modified", Covariate)) %>%
  mutate(Covariate=gsub("Primary.Lifestyle", "Primary lifestyle - ", Covariate)) %>%
  mutate(Covariate=gsub("Range.Size", "Range size", Covariate)) %>%
  mutate(Covariate=gsub("habitat_breadth", "Habitat breadth", Covariate)) %>%
  mutate(Covariate=gsub("Trophic.Level", "Trophic level - ", Covariate)) %>%
  mutate(Covariate=gsub("Mass", "Body mass", Covariate)) %>%
  arrange(Estimate) %>%
  ggplot(., aes(x=fct_inorder(Covariate), y=Estimate))+
  geom_point()+
  geom_errorbar(aes(x=fct_inorder(Covariate), y=Estimate, ymin=lwr, ymax=upr), width=0)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red", linetype="dashed")+
  coord_flip()+
  ylab("Parameter estimate")+
  xlab("")+
  ggtitle("Alpha Simpson")

ggsave("Figures/alpha_simpson_multiple_regression_results.png", width=6.9, height=6.4, units="in")


# plot the simpson vs richness values
mult_results_richness %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
  mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
  dplyr::select(Covariate, Estimate, upr, lwr) %>%
  rename(richness=Estimate,
         richness_upr=upr,
         richness_lwr=lwr) %>%
  left_join(., mult_results_simpson %>%
              mutate(Estimate=as.numeric(as.character(Estimate))) %>%
              mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
              mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
              dplyr::select(Covariate, Estimate, upr, lwr) %>%
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

ggsave("Figures/richness_vs_simpson_estimates_multiple_regression_results.png", height=5.6, width=6.1, units="in")


mult_results_richness %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
  mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
  dplyr::select(Covariate, Estimate, upr, lwr) %>%
  rename(richness=Estimate,
         richness_upr=upr,
         richness_lwr=lwr) %>%
  left_join(., mult_results_simpson %>%
              mutate(Estimate=as.numeric(as.character(Estimate))) %>%
              mutate(upr=as.numeric(as.character(`u-95% CI`))) %>%
              mutate(lwr=as.numeric(as.character(`l-95% CI`))) %>%
              dplyr::select(Covariate, Estimate, upr, lwr) %>%
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
  geom_smooth(method="lm")+
  geom_label_repel()









mult_results_richness_gamma <- readRDS("intermediate_results/gamma_multiple_regression_richness_summary.RDS") %>%
  mutate(scale="gamma") %>%
  mutate(response="richness")

mult_results_simpson_gamma <- readRDS("intermediate_results/gamma_multiple_regression_simpson_summary.RDS") %>%
  mutate(scale="gamma") %>%
  mutate(response="richness")
