# initial EDA
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(lme4)
library(ggeffects)
library(naniar)
library(performance)

# prepare trait data
analysis_dat <- readRDS("Clean data/analysis_data.RDS") %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

length(unique(analysis_dat$ebird_COMMON_NAME))
nrow(analysis_dat)
length(unique(analysis_dat$DOI))

# trim the data to traits of interest
# and the variables necessary for analyses
analysis_dat <- analysis_dat %>%
  dplyr::select(3, 7, 8, 9, 10, 12, 14, 15, 17, 20:22)

vis_miss(analysis_dat %>%
           dplyr::select(6:12))+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

analysis_dat %>%
  group_by(DOI) %>%
  summarize(number_obs=n(),
            number_species=length(unique(ebird_COMMON_NAME)))

summarized_dat <- analysis_dat %>%
  group_by(ebird_COMMON_NAME) %>%
  summarize(number_studies=length(unique(DOI)))

#############################
# body mass vs richness

# raw plots first
body_size_dat <- analysis_dat %>%
  dplyr::select(1, 3, 4, 6) %>%
  dplyr::filter(complete.cases(.)) %>%
  mutate(adult_body_mass_g_log10=log10(adult_body_mass_g)) %>%
  dplyr::filter(! ebird_COMMON_NAME %in% c("Emu", "Red Junglefowl"))

ggplot(body_size_dat, aes(x=adult_body_mass_g, y=Richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("log10 Adult body mass (g)")

ggplot(body_size_dat, aes(x=adult_body_mass_g, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("log10 Adult body mass (g)")

body_size_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + (1|DOI/ebird_COMMON_NAME),
                            family=poisson(), data=body_size_dat)
summary(body_size_mod)

body_size_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + (1|DOI) + (1|ebird_COMMON_NAME),
                             family=poisson(), data=body_size_dat)
summary(body_size_mod)

body_size_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + (1|DOI),
                             family=poisson(), data=body_size_dat)
summary(body_size_mod)

body_size_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + (1|ebird_COMMON_NAME),
                             family=poisson(), data=body_size_dat)
summary(body_size_mod)

hist(resid(body_size_mod))

model_performance(body_size_mod)

body_size_prediction <- ggpredict(body_size_mod, terms="adult_body_mass_g")

newdat <- data.frame(adult_body_mass_g=seq(min(body_size_dat$adult_body_mass_g_log10), max(body_size_dat$adult_body_mass_g_log10),
                                           length.out=100))

predict(body_size_mod, newdat, re.form=NA)

body_size_prediction <- data.frame(adult_body_mass_g_log10=newdat$adult_body_mass_g) %>%
  mutate(adult_body_mass_g=10^adult_body_mass_g_log10) %>%
  mutate(predicted_richness=exp(predict(body_size_mod, newdat, re.form=NA)))

ggplot(body_size_prediction, aes(x=adult_body_mass_g, y=predicted_richness)) +
  geom_line() +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=body_size_dat, aes(x=adult_body_mass_g, y=Richness))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Adult body mass")+
  scale_x_log10()

body_size_mod_l <- lme4::lmer(log10(Richness) ~ log10(adult_body_mass_g) + (1|DOI),
                             data=body_size_dat)
summary(body_size_mod_l)

body_size_prediction <- ggpredict(body_size_mod_l, terms="adult_body_mass_g")

ggplot(body_size_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  scale_x_log10()+
  geom_point(data=body_size_dat, aes(x=adult_body_mass_g, y=Richness))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Adult body mass")


body_size_mod_g <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + (1|ebird_COMMON_NAME),
                             family=poisson(), data=body_size_dat)
summary(body_size_mod_g)

model_performance(body_size_mod)

body_size_prediction <- ggpredict(body_size_mod_g, terms="adult_body_mass_g")

ggplot(body_size_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  scale_x_log10()+
  geom_point(data=body_size_dat, aes(x=adult_body_mass_g, y=Richness))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Adult body mass")

body_size_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + (1|ebird_COMMON_NAME),
                             family=poisson(), data=body_size_dat)
summary(body_size_mod)

model_performance(body_size_mod)

body_size_prediction <- ggpredict(body_size_mod, terms="adult_body_mass_g")

ggplot(body_size_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=body_size_dat, aes(x=adult_body_mass_g, y=Richness))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Adult body mass")

body_size_mod <- lme4::lmer(log10(Richness) ~ log10(adult_body_mass_g) + (1|DOI/ebird_COMMON_NAME),
                            data=body_size_dat)
summary(body_size_mod)

body_size_prediction <- ggpredict(body_size_mod, terms="adult_body_mass_g")

ggplot(body_size_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)

ggpredict(model = body_size_mod, terms="adult_body_mass_g", back.transform = TRUE) %>%
  plot(add.data = TRUE)

# get average richness per species?
body_size_dat_v2 <- body_size_dat %>%
  group_by(ebird_COMMON_NAME, DOI) %>%
  mutate(mean_richness=mean(Richness),
         sd_richness=sd(Richness)) %>%
  dplyr::select(-Richness) %>%
  distinct()

ggplot(body_size_dat_v2, aes(x=adult_body_mass_g, y=mean_richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("log10 Adult body mass (g)")

ggplot(body_size_dat_v2, aes(x=adult_body_mass_g, y=mean_richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("log10 Adult body mass (g)")+
  geom_smooth(method="lm")


#############################
# range size vs richness

# raw plots first
range_size_dat <- analysis_dat %>%
  dplyr::select(1, 3, 4, 12) %>%
  dplyr::filter(complete.cases(.))

ggplot(range_size_dat, aes(x=range_size, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("log10 Range size (km2)")

ggplot(range_size_dat, aes(x=range_size, y=Richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("log10 Range size (km2)")

hist(range_size_dat$Richness)

range_size_mod <- lme4::glmer(Richness ~ log10(range_size) + (1|DOI/ebird_COMMON_NAME),
                             family=poisson(), data=range_size_dat)
summary(range_size_mod)

range_size_prediction <- ggpredict(range_size_mod, terms="range_size")

ggplot(range_size_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)


#############################
# habitat breadth vs richness

# raw plots first
habitat_breadth_dat <- analysis_dat %>%
  dplyr::select(1, 3, 4, 11) %>%
  dplyr::filter(complete.cases(.))

hist(habitat_breadth_dat$habitat_breadth)

ggplot(habitat_breadth_dat, aes(x=habitat_breadth, y=Richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Habitat breadth")

ggplot(habitat_breadth_dat, aes(x=habitat_breadth, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("Habitat breadth (log10)")

habitat_breadth_mod <- lme4::glmer(Richness ~ habitat_breadth + (1|DOI/ebird_COMMON_NAME),
                              family=poisson(), data=habitat_breadth_dat)
summary(habitat_breadth_mod)

habitat_breadth_prediction <- ggpredict(habitat_breadth_mod, terms="habitat_breadth")

ggplot(habitat_breadth_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=habitat_breadth_dat, aes(x=habitat_breadth, y=Richness))+
  xlab("Habitat breadth")+
  ylab("Species richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

habitat_breadth_mod <- lme4::glmer(Richness ~ habitat_breadth + (1|DOI),
                                   family=poisson(), data=habitat_breadth_dat)
summary(habitat_breadth_mod)

habitat_breadth_prediction <- ggpredict(habitat_breadth_mod, terms="habitat_breadth")

ggplot(habitat_breadth_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=habitat_breadth_dat, aes(x=habitat_breadth, y=Richness))+
  xlab("Habitat breadth")+
  ylab("Species richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

habitat_breadth_mod <- lme4::glmer(Richness ~ habitat_breadth + (1|ebird_COMMON_NAME),
                                   family=poisson(), data=habitat_breadth_dat)
summary(habitat_breadth_mod)

habitat_breadth_prediction <- ggpredict(habitat_breadth_mod, terms="habitat_breadth")

ggplot(habitat_breadth_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=habitat_breadth_dat, aes(x=habitat_breadth, y=Richness))+
  xlab("Habitat breadth")+
  ylab("Species richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))




#############################
# flock size vs richness

# raw plots first
flock_size_dat <- analysis_dat %>%
  dplyr::select(1, 3, 4, 10) %>%
  dplyr::filter(complete.cases(.))

hist(flock_size_dat$mean_flock_size)
hist(log10(flock_size_dat$mean_flock_size))

ggplot(flock_size_dat, aes(x=mean_flock_size, y=Richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Mean flock size")

ggplot(flock_size_dat, aes(x=mean_flock_size, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("Mean flock size (log10)")

flock_size_mod <- lme4::glmer(Richness ~ mean_flock_size + (1|DOI/ebird_COMMON_NAME),
                                   family=poisson(), data=flock_size_dat)
summary(flock_size_mod)

flock_size_prediction <- ggpredict(flock_size_mod, terms="mean_flock_size")

ggplot(flock_size_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)


#############################
# diet breadth vs richness

# raw plots first
diet_breadth_dat <- analysis_dat %>%
  dplyr::select(1, 3, 4, 8) %>%
  dplyr::filter(complete.cases(.))

hist(diet_breadth_dat$diet_breadth)
hist(log10(diet_breadth_dat$diet_breadth))

ggplot(diet_breadth_dat, aes(x=diet_breadth, y=Richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness")+
  xlab("Diet breadth")

ggplot(diet_breadth_dat, aes(x=diet_breadth, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Species richness (log10)")+
  xlab("Diet breadth (log10)")

diet_breadth_mod <- lme4::glmer(Richness ~ diet_breadth + (1|DOI/ebird_COMMON_NAME),
                              family=poisson(), data=diet_breadth_dat)
summary(diet_breadth_mod)

diet_breadth_prediction <- ggpredict(diet_breadth_mod, terms="diet_breadth")

ggplot(diet_breadth_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=diet_breadth_dat, aes(x=diet_breadth, y=Richness))+
  xlab("Diet breadth")+
  ylab("Species richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

diet_breadth_mod <- lme4::glmer(Richness ~ diet_breadth + (1|DOI),
                                family=poisson(), data=diet_breadth_dat)
summary(diet_breadth_mod)

diet_breadth_prediction <- ggpredict(diet_breadth_mod, terms="diet_breadth")

ggplot(diet_breadth_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=diet_breadth_dat, aes(x=diet_breadth, y=Richness))+
  xlab("Diet breadth")+
  ylab("Species richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

diet_breadth_mod <- lme4::glmer(Richness ~ diet_breadth + (1|ebird_COMMON_NAME),
                                family=poisson(), data=diet_breadth_dat)
summary(diet_breadth_mod)

diet_breadth_prediction <- ggpredict(diet_breadth_mod, terms="diet_breadth")

ggplot(diet_breadth_prediction, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=diet_breadth_dat, aes(x=diet_breadth, y=Richness))+
  xlab("Diet breadth")+
  ylab("Species richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))




# All traits
full_traits <- analysis_dat %>%
  dplyr::select(1:6, 8, 10:12) %>%
  dplyr::filter(complete.cases(.))

big_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + diet_breadth + habitat_breadth + mean_flock_size +
                         (1|DOI/ebird_COMMON_NAME), family=poisson(), data=full_traits)

summary(big_mod)

broom.mixed::tidy(big_mod) %>%
  slice(1:5) %>%
  mutate(lwr_95=confint(big_mod) %>%
           as.data.frame() %>%
           slice(2:6) %>%
           .[,1]) %>%
  mutate(upr_95=confint(big_mod) %>%
           as.data.frame() %>%
           slice(2:6) %>%
           .[,2]) %>%
  ggplot(., aes(x=term, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lwr_95, ymax=upr_95), width=0)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

big_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + diet_breadth + habitat_breadth + mean_flock_size +
                         (1|DOI), family=poisson(), data=full_traits)

broom.mixed::tidy(big_mod) %>%
  slice(1:5) %>%
  mutate(lwr_95=confint(big_mod) %>%
           as.data.frame() %>%
           slice(2:6) %>%
           .[,1]) %>%
  mutate(upr_95=confint(big_mod) %>%
           as.data.frame() %>%
           slice(2:6) %>%
           .[,2]) %>%
  ggplot(., aes(x=term, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lwr_95, ymax=upr_95), width=0)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

big_mod <- lme4::glmer(Richness ~ log10(adult_body_mass_g) + diet_breadth + habitat_breadth + mean_flock_size +
                         (1|ebird_COMMON_NAME), family=poisson(), data=full_traits)

broom.mixed::tidy(big_mod) %>%
  slice(1:5) %>%
  mutate(lwr_95=confint(big_mod) %>%
           as.data.frame() %>%
           slice(2:6) %>%
           .[,1]) %>%
  mutate(upr_95=confint(big_mod) %>%
           as.data.frame() %>%
           slice(2:6) %>%
           .[,2]) %>%
  ggplot(., aes(x=term, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lwr_95, ymax=upr_95), width=0)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  coord_flip()+
  geom_hline(yintercept=0, color="red", linetype="dashed")

summary(big_mod)

plot(big_mod)
hist(resid(big_mod))











ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=adult_body_mass_g, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=adult_body_mass_g, y=Richness))+
  geom_point()+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=adult_body_mass_g, y=ISimpson))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=adult_body_mass_g, y=Shannon))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=brain_residual, y=Richness))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=diet_breadth, y=Richness))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=diet_breadth, y=ISimpson))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=TrophicLevel, y=ISimpson))+
  geom_violin()+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=TrophicLevel, y=Richness))+
  geom_violin()+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=TrophicLevel, y=Shannon))+
  geom_violin()+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=TrophicNiche, y=ISimpson))+
  geom_violin()+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=TrophicNiche, y=Richness))+
  geom_violin()+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat %>%
         dplyr::filter(!ebird_COMMON_NAME %in% c("Red Junglefowl", "Emu")),
       aes(x=TrophicNiche, y=Shannon))+
  geom_violin()+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat, aes(x=Habitat_Breadth, y=Richness))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")

ggplot(analysis_dat, aes(x=Habitat_Breadth, y=ISimpson))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")
