## concatenate data sources by species
library(dplyr)
library(taxonlookup)

# read in data files
leafmorph <- read_csv("data/netobradley_garden_data.csv")
leafmetab <- read_csv("data/ecophysfit_thinned15.csv")
leafchemi <- read_csv("data/neto-bradley_msc_leaf_chemistry.csv")

# clean up files by species
leafchemi <- leafchemi %>% select(1:4,)

leafmetab$binomial <- gsub('^.{16}', '', leafmetab$UniqueID)
leafmetab$binomial <- str_to_title(leafmetab$binomial)

#create new dataframes with mean and sd by species
leafmetab_clean <- leafmetab %>% group_by(binomial) %>% summarise(vcmax_mean = mean(Vcmax), vcmax_sd = sd(Vcmax), jmax_mean = mean(Jmax), jmax_sd = sd(Jmax))
number <- leafmetab %>% group_by(binomial) %>% tally()
leafmetab_clean <- inner_join(leafmetab_clean, number)

leafmorph$thickmean <- (leafmorph$thick_1+leafmorph$thick_2+leafmorph$thick_3)/3
leafmorph_clean <- leafmorph %>% group_by(binomial) %>% summarise(thick_mean = mean(thickmean), thick_sd = sd(thickmean), area_mean = mean(leaf_area_imagej), log_area_mean = mean(log(leaf_area_imagej)), area_sd = sd(leaf_area_imagej), log_area_sd = sd(log(leaf_area_imagej)), freshmass_mean = mean(fresh_mass), freshmass_sd = sd(fresh_mass), drymass_mean = mean(dry_mass), drymass_sd = sd(dry_mass), log_freshmass_mean = mean(log(fresh_mass)), log_freshmass_sd = sd(log(fresh_mass)), sla_mean = mean(log(leaf_area_imagej/fresh_mass)), sla_sd = sd(log(leaf_area_imagej/fresh_mass)))

## remove Rhododendron campylocarpum
which(leafmetab_clean$binomial == "Rhododendron_campylocarpum")
#leafmetab_clean <- leafmetab_clean[-101,]

#write_csv(leafchemi, "data/cleanleafchemi.csv")
#write_csv(leafmetab_clean, "data/cleanleafmetab_thinned15.csv")
#write_csv(leafmorph_clean, "data/cleanleafmorph.csv")
#write_csv(leafmorph_clean, "data/cleanleafmorphnormal.csv")
