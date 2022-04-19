# a script to make figures of the modelled results

# packages
library(dplyr)
library(ggplot2)
library(ggstance)
library(tidybayes)
library(tidyr)

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
                                      "mean_flock_size", "Habitat", "Trophic.Level", "Trophic.Niche",
                                      "Primary.Lifestyle"), get_summary))

############################
############################
# make a figure for the Mass
############################
############################
mass_draws <- readRDS("intermediate_results/single_regression_richness_draws_Mass.RDS")
mass_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_Mass.RDS")
mass_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_Mass.RDS")
mass_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, Mass) %>%
  dplyr::filter(complete.cases(.))

mass_plot_raw <- ggplot()+
  geom_point(data=mass_data, aes(x=Mass, y=Richness), alpha=0.4)+
  geom_line(data=mass_fe_only, aes(x=predictor, y=.value, group=.draw), 
            alpha=0.1)+
  geom_line(data=mass_fe_only_mean, aes(x=predictor, y=.value),
            color="red", lwd=2, group=1)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Body mass (g)")+
  ylab("Microbial species richness")+
  ggtitle("A")

mass_plot_raw

mass_plot_posterior <- ggplot()+
  stat_halfeye(data=mass_draws, aes(x=b_log10predictor))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  xlab("Posterior distribution")+
  ylab("")

mass_plot_posterior

##################################
##################################
# make a figure for the Range size
##################################
##################################
Range.Size_draws <- readRDS("intermediate_results/single_regression_richness_draws_Range.Size.RDS")
Range.Size_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_Range.Size.RDS")
Range.Size_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_Range.Size.RDS")
Range.Size_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, Range.Size) %>%
  dplyr::filter(complete.cases(.))

Range.Size_plot <- ggplot()+
  geom_point(data=Range.Size_data, aes(x=Range.Size, y=Richness), alpha=0.4)+
  geom_line(data=Range.Size_fe_only, aes(x=predictor, y=.value, group=.draw), 
            alpha=0.05)+
  geom_line(data=Range.Size_fe_only_mean, aes(x=predictor, y=.value),
            color="red", lwd=2, group=1)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Range size (km2)")+
  ylab("Microbial species richness")+
  ggtitle("B")

Range.Size_plot

Range.Size_plot_posterior <- ggplot()+
  stat_halfeye(data=Range.Size_draws, aes(x=b_log10predictor))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  xlab("Posterior distribution")+
  ylab("")

Range.Size_plot_posterior

########################################
########################################
# make a figure for the habitat breadth
########################################
########################################
habitat_breadth_draws <- readRDS("intermediate_results/single_regression_richness_draws_habitat_breadth.RDS")
habitat_breadth_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_habitat_breadth.RDS")
habitat_breadth_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_habitat_breadth.RDS")
habitat_breadth_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, habitat_breadth) %>%
  dplyr::filter(complete.cases(.))

habitat_breadth_plot <- ggplot()+
  geom_point(data=habitat_breadth_data, aes(x=habitat_breadth, y=Richness), alpha=0.4)+
  geom_line(data=habitat_breadth_fe_only, aes(x=predictor, y=.value, group=.draw), 
            alpha=0.05)+
  geom_line(data=habitat_breadth_fe_only_mean, aes(x=predictor, y=.value),
            color="red", lwd=2, group=1)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Habitat breadth")+
  ylab("Microbial species richness")+
  ggtitle("C")

habitat_breadth_plot

habitat_breadth_plot_posterior <- ggplot()+
  stat_halfeye(data=habitat_breadth_draws, aes(x=b_log10predictor))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  xlab("Posterior distribution")+
  ylab("")

habitat_breadth_plot_posterior

###################################
###################################
# make a figure for the Flock size
###################################
###################################
flock_size_draws <- readRDS("intermediate_results/single_regression_richness_draws_mean_flock_size.RDS")
flock_size_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_flock_size.RDS")
flock_size_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_mean_flock_size.RDS")
flock_size_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, mean_flock_size) %>%
  dplyr::filter(complete.cases(.))

flock_size_plot_raw <- ggplot()+
  geom_point(data=flock_size_data, aes(x=mean_flock_size, y=Richness), alpha=0.4)+
  geom_line(data=flock_size_fe_only, aes(x=predictor, y=.value, group=.draw), 
            alpha=0.1)+
  geom_line(data=flock_size_fe_only_mean, aes(x=predictor, y=.value),
            color="red", lwd=2, group=1)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Mean flock size)")+
  ylab("Microbial species richness")+
  ggtitle("A")

flock_size_plot_raw

flock_size_plot_posterior <- ggplot()+
  stat_halfeye(data=flock_size_draws, aes(x=b_log10predictor))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  xlab("Posterior distribution")+
  ylab("")

flock_size_plot_posterior

###################################
###################################
# make a figure for the Abundance
###################################
###################################
pop_abund_draws <- readRDS("intermediate_results/single_regression_richness_draws_pop_abund.RDS")
pop_abund_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_pop_abund.RDS")
pop_abund_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_pop_abund.RDS")
pop_abund_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, pop_abund) %>%
  dplyr::filter(complete.cases(.))

pop_abund_plot_raw <- ggplot()+
  geom_point(data=pop_abund_data, aes(x=pop_abund, y=Richness), alpha=0.4)+
  geom_line(data=pop_abund_fe_only, aes(x=predictor, y=.value, group=.draw), 
            alpha=0.1)+
  geom_line(data=pop_abund_fe_only_mean, aes(x=predictor, y=.value),
            color="red", lwd=2, group=1)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Population abundance")+
  ylab("Microbial species richness")+
  ggtitle("A")

pop_abund_plot_raw

pop_abund_plot_posterior <- ggplot()+
  stat_halfeye(data=pop_abund_draws, aes(x=b_log10predictor))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  xlab("Posterior distribution")+
  ylab("")

pop_abund_plot_posterior

###################################
###################################
# make a figure for the Habitat type
###################################
###################################
habitat_draws <- readRDS("intermediate_results/single_regression_richness_draws_Habitat.RDS")
habitat_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_Habitat.RDS")
habitat_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_Habitat.RDS")
habitat_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, Habitat) %>%
  dplyr::filter(complete.cases(.))

habitat_plot_raw <- habitat_fe_only %>%
  ggplot(aes(x=predictor, y=.value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Habitat type")+
  ylab("Microbial species richness")+
  ggtitle("A")+
  coord_flip()

habitat_plot_raw

habitat_plot_posterior <- habitat_draws %>%
  pivot_longer(!c(1:3), names_to="term", values_to="value") %>%
  ggplot(data=., aes(x=term, y=value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ylab("Posterior distribution")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

habitat_plot_posterior

###################################
###################################
# make a figure for the trophic niche
###################################
###################################
trophic_niche_draws <- readRDS("intermediate_results/single_regression_richness_draws_Trophic.Niche.RDS")
trophic_niche_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_Trophic.Niche.RDS")
trophic_niche_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_Trophic.Niche.RDS")
trophic_niche_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, Trophic.Niche) %>%
  dplyr::filter(complete.cases(.))

trophic_niche_plot_raw <- trophic_niche_fe_only %>%
  ggplot(aes(x=predictor, y=.value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Trophic niche")+
  ylab("Microbial species richness")+
  ggtitle("A")+
  coord_flip()

trophic_niche_plot_raw

trophic_niche_plot_posterior <- trophic_niche_draws %>%
  pivot_longer(!c(1:3), names_to="term", values_to="value") %>%
  ggplot(data=., aes(x=term, y=value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ylab("Posterior distribution")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

trophic_niche_plot_posterior

###################################
###################################
# make a figure for the trophic level
###################################
###################################
trophic_level_draws <- readRDS("intermediate_results/single_regression_richness_draws_Trophic.Level.RDS")
trophic_level_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_Trophic.Level.RDS")
trophic_level_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_Trophic.Level.RDS")
trophic_level_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, Trophic.Level) %>%
  dplyr::filter(complete.cases(.))

trophic_level_plot_raw <- trophic_level_fe_only %>%
  ggplot(aes(x=predictor, y=.value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Trophic niche")+
  ylab("Microbial species richness")+
  ggtitle("A")+
  coord_flip()

trophic_level_plot_raw

trophic_level_plot_posterior <- trophic_level_draws %>%
  pivot_longer(!c(1:3), names_to="term", values_to="value") %>%
  ggplot(data=., aes(x=term, y=value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ylab("Posterior distribution")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

trophic_level_plot_posterior

###################################
###################################
# make a figure for the primary lifestyle
###################################
###################################
primary_lifestyle_draws <- readRDS("intermediate_results/single_regression_richness_draws_Primary.Lifestyle.RDS")
primary_lifestyle_fe_only <- readRDS("intermediate_results/single_regression_richness_fe_only_Primary.Lifestyle.RDS")
primary_lifestyle_fe_only_mean <- readRDS("intermediate_results/single_regression_richness_fe_only_mean_Primary.Lifestyle.RDS")
primary_lifestyle_data <- analysis_dat %>%
  dplyr::select(ebird_COMMON_NAME, Richness, ISimpson, DOI, Sample_type, Primary.Lifestyle) %>%
  dplyr::filter(complete.cases(.))

primary_lifestyle_plot_raw <- primary_lifestyle_fe_only %>%
  ggplot(aes(x=predictor, y=.value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Primary lifestyle")+
  ylab("Microbial species richness")+
  ggtitle("A")+
  coord_flip()

primary_lifestyle_plot_raw

primary_lifestyle_plot_posterior <- primary_lifestyle_draws %>%
  pivot_longer(!c(1:3), names_to="term", values_to="value") %>%
  ggplot(data=., aes(x=term, y=value))+
  stat_halfeye()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ylab("Posterior distribution")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

primary_lifestyle_plot_posterior

