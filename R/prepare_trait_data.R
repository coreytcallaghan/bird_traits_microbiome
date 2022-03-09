# this script is used to prepare the 'predictor' variables (i.e., life history traits)
# in order to get them ready for analysis
# it uses the 'metadata' file from Steph
# but the 'diversity' file does differ slightly from this metadata file
# nevertheless, this is cleaned a little bit further in the
# prepare_data_for_analysis.R script

# packages
library(dplyr)
library(tidyr)
library(readr)
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

# now do it again
species_list_clements2 <- species %>%
  left_join(., fixes) %>%
  mutate(ebird_SCIENTIFIC_NAME=ifelse(is.na(ebird_SCIENTIFIC_NAME)==TRUE, 
                                      SCIENTIFIC_NAME, ebird_SCIENTIFIC_NAME)) %>%
  dplyr::select(SCIENTIFIC_NAME, ebird_SCIENTIFIC_NAME) %>%
  left_join(., clements, by="ebird_SCIENTIFIC_NAME") %>%
  rename(database_SCIENTIFIC_NAME=SCIENTIFIC_NAME)

# now start joining trait data
# BODY SIZE 
body_size <- read_csv("trait_data/body size data/cleaned_body_size_data.csv") %>%
  rename(ebird_COMMON_NAME=COMMON_NAME) %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME)

body_dat1 <- species_list_clements2 %>%
  left_join(., body_size %>%
              dplyr::select(1, 7)) %>%
  dplyr::filter(complete.cases(adult_body_mass_g))

body_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% body_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., body_size %>%
              dplyr::select(2, 7))

body_dat_full <- body_dat1 %>%
  bind_rows(body_dat2) %>%
  dplyr::filter(complete.cases(adult_body_mass_g))

body_dat3 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% body_dat_full$ebird_SCIENTIFIC_NAME) %>%
  left_join(., body_size %>%
              dplyr::select(3, 7))

body_dat_full <- body_dat_full %>%
  bind_rows(body_dat3)

# BRAIN SIZE 
brain_size <- read_csv("trait_data/brain_size_data/brain_size_and_other_data.csv") %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME) %>%
  mutate(TipLabel=ebird_SCIENTIFIC_NAME) %>%
  mutate(ebird_SCIENTIFIC_NAME=gsub("_", " ", ebird_SCIENTIFIC_NAME))

brain_dat1 <- species_list_clements2 %>%
  left_join(., brain_size %>%
              dplyr::select(1, 4:10)) %>%
  dplyr::filter(complete.cases(brain_size_mm3))

brain_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% brain_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., brain_size %>%
              dplyr::select(4:11))

brain_dat_full <- brain_dat1 %>%
  bind_rows(brain_dat2)

# ecological niche assignment
niche <- read_csv("trait_data/ecological_niche_assignment/41559_2019_1070_MOESM3_ESM.csv") %>%
  rename(ebird_SCIENTIFIC_NAME=Binomial) %>%
  mutate(TipLabel=ebird_SCIENTIFIC_NAME) %>%
  mutate(ebird_SCIENTIFIC_NAME=gsub("_", " ", ebird_SCIENTIFIC_NAME))

niche_dat1 <- species_list_clements2 %>%
  left_join(., niche %>%
              dplyr::select(1,15:18)) %>%
  dplyr::filter(complete.cases(TrophicNiche))

niche_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% niche_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., niche %>%
              dplyr::select(15:18, 20))

niche_dat_full <- niche_dat1 %>%
  bind_rows(niche_dat2)

# flock size (gregariousness)
flock_size <- readRDS("trait_data/flock size/flock_size_per_month_for_each_species.RDS") %>%
  group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(mean_flock_size=mean(mean_abund)) %>%
  mutate(ebird_COMMON_NAME=COMMON_NAME) %>%
  mutate(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME)

flock_dat1 <- species_list_clements2 %>%
  left_join(., flock_size %>%
              dplyr::select(3, 5)) %>%
  dplyr::filter(complete.cases(mean_flock_size))

flock_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% flock_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., flock_size %>%
              dplyr::select(3, 4))

flock_dat_full <- flock_dat1 %>%
  bind_rows(flock_dat2) %>%
  dplyr::select(-COMMON_NAME)

# fecundity and life history
clutch_size <- read_delim("trait_data/fecundity_and_life_history/avian_ssd_jan07.txt", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(ebird_COMMON_NAME=English_name) %>%
  mutate(ebird_SCIENTIFIC_NAME=Species_name)

clutch_dat1 <- species_list_clements2 %>%
  left_join(., clutch_size %>%
              dplyr::select(36, 46)) %>%
  dplyr::filter(complete.cases(Clutch_size))

clutch_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% clutch_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., clutch_size %>%
              dplyr::select(36, 45))

clutch_dat_full <- clutch_dat1 %>%
  bind_rows(clutch_dat2) %>%
  dplyr::filter(complete.cases(Clutch_size))

clutch_dat3 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% clutch_dat_full$ebird_SCIENTIFIC_NAME) %>%
  left_join(., clutch_size %>%
              mutate(TipLabel=gsub(" ", "_", ebird_SCIENTIFIC_NAME)) %>%
              dplyr::select(36, 47))

clutch_dat_full <- clutch_dat_full %>%
  bind_rows(clutch_dat3) %>%
  mutate(Clutch_size=ifelse(Clutch_size==-999, NA, Clutch_size)) %>%
  distinct()

# migration status
migration <- read_csv("trait_data/migration_status/migration_status_cleaned.csv") %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, MIGRATORY_STATUS) %>%
  distinct() %>%
  mutate(ebird_COMMON_NAME=COMMON_NAME) %>%
  mutate(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME)

mig_dat1 <- species_list_clements2 %>%
  left_join(., migration %>%
              dplyr::select(3, 5)) %>%
  dplyr::filter(complete.cases(MIGRATORY_STATUS))

mig_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% mig_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., migration %>%
              dplyr::select(3, 4))

mig_dat_full <- mig_dat1 %>%
  bind_rows(mig_dat2)

# iucn habitats
iucn <- read_csv("trait_data/iucn_habitats/iucn_habitats.csv") %>%
  dplyr::select(8, 68) %>%
  rename(ebird_SCIENTIFIC_NAME=Sciname) %>%
  mutate(TipLabel=gsub(" ", "_", ebird_SCIENTIFIC_NAME))

iucn_dat1 <- species_list_clements2 %>%
  left_join(., iucn %>%
              dplyr::select(1, 2)) %>%
  dplyr::filter(complete.cases(Habitat_Breadth))

iucn_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% iucn_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., iucn %>%
              dplyr::select(2, 3))

iucn_dat_full <- iucn_dat1 %>%
  bind_rows(iucn_dat2) %>%
  dplyr::filter(complete.cases(Habitat_Breadth))

# range size
range <-  read_delim("trait_data/range_size/bird.range.season041219.csv", 
                     ";", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(1, 6) %>%
  rename(ebird_SCIENTIFIC_NAME=Species) %>%
  mutate(TipLabel=gsub(" ", "_", ebird_SCIENTIFIC_NAME))

range_dat1 <- species_list_clements2 %>%
  left_join(., range %>%
              dplyr::select(1, 2)) %>%
  dplyr::filter(complete.cases(range.size.km2))

range_dat2 <- species_list_clements2 %>%
  dplyr::filter(! ebird_SCIENTIFIC_NAME %in% range_dat1$ebird_SCIENTIFIC_NAME) %>%
  left_join(., range %>%
              dplyr::select(2, 3))

range_dat_full <- range_dat1 %>%
  bind_rows(range_dat2) %>%
  dplyr::filter(complete.cases(range.size.km2))

# avonet traits
# read in trait data
avonet_traits <- read_csv("AVONET_data/Supplementary dataset 1.csv") %>%
  rename(ebird_SCIENTIFIC_NAME=Species2)

avonet_dat <- species_list_clements2 %>%
  left_join(., avonet_traits)

# combine into one dataframe
trait_dat <- body_dat_full %>%
  left_join(., brain_dat_full) %>%
  left_join(., niche_dat_full) %>%
  left_join(., flock_dat_full) %>%
  left_join(., clutch_dat_full) %>%
  left_join(., mig_dat_full) %>%
  left_join(., iucn_dat_full) %>%
  left_join(., range_dat_full) %>%
  left_join(., avonet_dat)

saveRDS(trait_dat, "Clean data/bird_trait_predictors.RDS")




