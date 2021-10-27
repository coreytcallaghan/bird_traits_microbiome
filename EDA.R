# initial EDA

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)

# prepare trait data
dat <- readr::read_csv("All_Metadata.csv")

# get unique list of species
species <- dat %>%
  mutate(Species2=word(Species, 2, sep="_")) %>%
  mutate(Species3=ifelse(is.na(Species2)==TRUE, Species, Species2)) %>%
  unite(SCIENTIFIC_NAME, Genus, Species3, sep=" ") %>%
  dplyr::filter(complete.cases(Species)) %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  dplyr::filter(SCIENTIFIC_NAME != "Anser sp.")

clements <- read_csv("clements_clean.csv")

species_list_clements <- species %>%
  left_join(., clements %>%
              rename(SCIENTIFIC_NAME=ebird_SCIENTIFIC_NAME))

fixes <- species_list_clements %>%
  dplyr::filter(is.na(ebird_COMMON_NAME)) %>%
  mutate(ebird_SCIENTIFIC_NAME=c("Paramythia montium",
                                 "Melipotes fumigatus",
                                 "NA",
                                 "Arses insularis",
                                 "Symposiachrus guttula",
                                 "Turtur tympanistria",
                                 "Corythornis leucogaster",
                                 "Cyanomitra olivacea",
                                 "Spermestes bicolor",
                                 "Spermestes cucullata",
                                 "Cecropis abyssinica",
                                 "Alethe castanea",
                                 "Chamaetylas poliocephala",
                                 "Bradornis fuliginosus",
                                 "Hedydipna collaris",
                                 "Anthreptes seimundi",
                                 "Cinnyris batesi",
                                 "Cinnyris reichenowi",
                                 "Cyanomitra oritis",
                                 "Arizelocichla tephrolaema",
                                 "Eurillas latirostris",
                                 "Eurillas virens",
                                 "Sylvia abyssinica",
                                 "Zosterops brunneus",
                                 "Gallus gallus",
                                 "Geospiza difficilis",
                                 "Geospiza fortis",
                                 "Geospiza fuliginosa",
                                 "Geospiza scandens",
                                 "Camarhynchus parvulus",
                                 "Platyspiza crassirostris",
                                 "Geospiza conirostris",
                                 "Certhidea olivacea",
                                 "Camarhynchus pallidus",
                                 "Turdus flavipes"))


alpha <- readr::read_csv("AlphaDiversity.csv")

# prepare trait data
dat <- readr::read_csv("All_Metadata.csv") %>%
  dplyr::filter(SampleID %in% alpha$SampleID)

predictors <- readRDS("bird_trait_predictors.RDS")


# get unique list of species
species <- alpha %>%
  mutate(Species2=word(Species, 2, sep="_")) %>%
  mutate(Species3=ifelse(is.na(Species2)==TRUE, Species, Species2)) %>%
  unite(SCIENTIFIC_NAME, Genus, Species3, sep=" ") %>%
  dplyr::filter(complete.cases(Species)) %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  dplyr::filter(SCIENTIFIC_NAME != "Anser sp.")

alpha2 <- species %>%
  left_join(., fixes) %>%
  mutate(ebird_SCIENTIFIC_NAME=ifelse(is.na(ebird_SCIENTIFIC_NAME)==TRUE, 
                                      SCIENTIFIC_NAME, ebird_SCIENTIFIC_NAME)) %>%
  dplyr::select(SCIENTIFIC_NAME, ebird_SCIENTIFIC_NAME) %>%
  left_join(., clements, by="ebird_SCIENTIFIC_NAME") %>%
  rename(database_SCIENTIFIC_NAME=SCIENTIFIC_NAME) %>%
  left_join(., alpha %>%
              mutate(Species2=word(Species, 2, sep="_")) %>%
              mutate(Species3=ifelse(is.na(Species2)==TRUE, Species, Species2)) %>%
              unite(database_SCIENTIFIC_NAME, Genus, Species3, sep=" "))

analysis_dat <- alpha2 %>%
  left_join(., predictors)


ggplot(analysis_dat, aes(x=adult_body_mass_g, y=Richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")+
  geom_smooth(method="gam")

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
