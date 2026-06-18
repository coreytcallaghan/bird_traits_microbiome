# a script to see what to do for the new analysis
# going to try and understand the new files Steph sent here.

library(dplyr)
library(readr)

# compare old and new files for richness
gamma_old <- read.table("Raw data/gamma_Richness.txt", header = TRUE)
gamma_new_raw <- read.table("Corey2026/Gammaq0.txt", header = TRUE) %>%
  mutate(NameXSpecies=row.names(.))

# Compare row identifiers
old_ids <- gamma_old$NameXSpecies
new_ids <- gamma_new_raw$NameXSpecies

missing_from_new <- setdiff(old_ids, new_ids)
new_not_in_old <- setdiff(new_ids, old_ids)

missing_from_new
new_not_in_old

length(missing_from_new)
length(new_not_in_old)

# Look at full old rows that disappeared
gamma_old %>%
  filter(NameXSpecies %in% missing_from_new)


gamma_richness_v2 <- gamma_new_raw %>%
  rowwise() %>%
  mutate(
    mean = mean(c_across(starts_with("V")), na.rm = TRUE),
    sd = sd(c_across(starts_with("V")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    ScientificName = str_extract(NameXSpecies, "^[A-Za-z]+\\s+[^ ]+"),
    Study = str_extract(NameXSpecies, "10\\.[^ ]+|Missing: Not provided"),
    Site = str_trim(str_remove(NameXSpecies, paste0("^", ScientificName, "\\s+", Study, "\\s+")))
  ) %>%
  select(NameXSpecies, mean, sd, ScientificName, Study, Site)


# now do the same as above but for ISimpson
gamma_old_simpson <- read.table("Raw data/gamma_ISimpson.txt", header = TRUE)
gamma_new_raw_simpson <- read.table("Corey2026/Gammaq2.txt", header = TRUE) %>%
  mutate(NameXSpecies=row.names(.))

# Compare row identifiers
old_ids_simpson <- gamma_old_simpson$NameXSpecies
new_ids_simpson <- gamma_new_raw_simpson$NameXSpecies

missing_from_new_simpson <- setdiff(old_ids_simpson, new_ids_simpson)
new_not_in_old_simpson <- setdiff(new_ids_simpson, old_ids_simpson)

missing_from_new_simpson
new_not_in_old_simpson

length(missing_from_new_simpson)
length(new_not_in_old_simpson)

# Look at full old rows that disappeared
gamma_old_simpson %>%
  filter(NameXSpecies %in% missing_from_new_simpson)


gamma_simpson_v2 <- gamma_new_raw_simpson %>%
  rowwise() %>%
  mutate(
    mean = mean(c_across(starts_with("V")), na.rm = TRUE),
    sd = sd(c_across(starts_with("V")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    ScientificName = str_extract(NameXSpecies, "^[A-Za-z]+\\s+[^ ]+"),
    Study = str_extract(NameXSpecies, "10\\.[^ ]+|Missing: Not provided"),
    Site = str_trim(str_remove(NameXSpecies, paste0("^", ScientificName, "\\s+", Study, "\\s+")))
  ) %>%
  select(NameXSpecies, mean, sd, ScientificName, Study, Site)



# Now that we've done the above, we pull old code from "prepare_data_for_analysis_gamma.R" and put it here
# and adapt it slightly to get to the same point (since we already have the gamma files above now)

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
analysis_dat_cleaned_richness <- analysis_dat %>%
  left_join(., gamma_richness_v2 %>%
              rename(database_SCIENTIFIC_NAME=ScientificName) %>%
              rename(DOI=Study) %>%
              rename(Sample_type=Site) %>%
              dplyr::select(2:4)) %>%
  dplyr::select(1:6, 20, 22, 42, 45, 65, 71:76, 84, 85, 86, 87) %>%
  distinct() %>%
  rename(habitat_breadth=Habitat_Breadth) %>%
  mutate(Migration=case_when(Migration==1 ~ "Sedentary",
                             Migration==2 ~ "Partially migratory",
                             Migration==3 ~ "Migratory")) %>%
  dplyr::filter(complete.cases(Mass)) %>%
  mutate(ebird_COMMON_NAME=ifelse(ebird_SCIENTIFIC_NAME=="Geospiza conirostris", "Espanola Ground-Finch", ebird_COMMON_NAME)) %>%
  dplyr::filter(complete.cases(mean)) %>%
  rename(gamma_Richness=mean,
         gamma_Richness_sd=sd)

analysis_dat_cleaned_isimpson <- analysis_dat %>%
  left_join(., gamma_simpson_v2 %>%
              rename(database_SCIENTIFIC_NAME=ScientificName) %>%
              rename(DOI=Study) %>%
              rename(Sample_type=Site) %>%
              dplyr::select(2:4)) %>%
  dplyr::select(1:6, 20, 22, 42, 45, 65, 71:76, 84, 85, 86, 87) %>%
  distinct() %>%
  rename(habitat_breadth=Habitat_Breadth) %>%
  mutate(Migration=case_when(Migration==1 ~ "Sedentary",
                             Migration==2 ~ "Partially migratory",
                             Migration==3 ~ "Migratory")) %>%
  dplyr::filter(complete.cases(Mass)) %>%
  mutate(ebird_COMMON_NAME=ifelse(ebird_SCIENTIFIC_NAME=="Geospiza conirostris", "Espanola Ground-Finch", ebird_COMMON_NAME)) %>%
  dplyr::filter(complete.cases(mean)) %>%
  rename(gamma_ISimpson=mean,
         gamma_ISimpson_sd=sd)

analysis_dat_cleaned <- analysis_dat_cleaned_richness %>%
  left_join(., analysis_dat_cleaned_isimpson)

saveRDS(analysis_dat_cleaned, "Clean data/analysis_data_gamma_v2.RDS")
