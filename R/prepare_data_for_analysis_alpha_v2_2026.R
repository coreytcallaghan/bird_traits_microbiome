# This script is intended to
# pull all the data together
# and create one output 'analysis' dataframe
# that will be used for visualization and
# analyses
# it relies on the 'predictor_variables.RDS' file outputted from the "prepare_trait_data.R" script
# but also because the alpha diversity file changes from the metadata originally sent
# there is some more hacking and manual fixes of taxonomy

# packages
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)

# prepare trait data
dat <- readr::read_csv("Raw data/All_Metadata.csv")

# get unique list of species
species <- dat %>%
  mutate(Species2=word(Species, 2, sep="_")) %>%
  mutate(Species3=ifelse(is.na(Species2)==TRUE, Species, Species2)) %>%
  unite(SCIENTIFIC_NAME, Genus, Species3, sep=" ") %>%
  dplyr::filter(complete.cases(Species)) %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  dplyr::filter(SCIENTIFIC_NAME != "Anser sp.")

clements <- read_csv("Raw data/clements_clean.csv")

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


alpha <- readr::read_csv("Raw data/AlphaDiversity.csv")


# remove old data no longer used.
removed_gamma_ids <- c(
  "Eurillas (Andropadus) virens 10.1111/mec.15354 intestine",
  "Paramythia montium 10.1098/rspb.2021.0446 Cloacal swabs",
  "Sylvia (Pseudoalcippe) abyssinica 10.1111/mec.15354 intestine"
)

alpha_check <- alpha %>%
  mutate(
    Species2 = word(Species, 2, sep = "_"),
    Species3 = ifelse(is.na(Species2), Species, Species2),
    NameXSpecies = paste(Genus, Species3, DOI, Sample_type)
  )

# they do exist
alpha_check %>%
  filter(NameXSpecies %in% removed_gamma_ids) %>%
  count(NameXSpecies)

# so now remove them
alpha_v2 <- alpha_check %>%
  filter(!NameXSpecies %in% removed_gamma_ids) %>%
  select(-Species2, -Species3, -NameXSpecies)

# rename alpha alpha v2 (the new, smaller file)
alpha <- alpha_v2

# then everything from here on stays the same
# prepare trait data
dat <- readr::read_csv("Raw data/All_Metadata.csv") %>%
  dplyr::filter(SampleID %in% alpha$SampleID)

predictors <- readRDS("Clean data/bird_trait_predictors.RDS")


# get unique list of species
# again
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

# now clean up the data
# and drop some columns
analysis_dat_cleaned <- analysis_dat %>%
  dplyr::filter(complete.cases(Richness)) %>%
  dplyr::select(1:6, 16, 20, 22, 25:27, 42, 45, 65, 71:76, 84, 85) %>%
  rename(habitat_breadth=Habitat_Breadth) %>%
  mutate(Migration=case_when(Migration==1 ~ "Sedentary",
                             Migration==2 ~ "Partially migratory",
                             Migration==3 ~ "Migratory")) %>%
  dplyr::filter(complete.cases(Mass)) %>%
  mutate(ebird_COMMON_NAME=ifelse(ebird_SCIENTIFIC_NAME=="Geospiza conirostris", "Espanola Ground-Finch", ebird_COMMON_NAME))

saveRDS(analysis_dat_cleaned, "Clean data/analysis_data_alpha_v2.RDS")
