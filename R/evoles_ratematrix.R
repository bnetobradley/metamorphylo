## packages used
library(ape)
library(devtools)
library(dplyr)
library(geiger)
library(phangorn)
library(phytools)
library(ratematrix)
library(taxonlookup)

## data used
all_plants_tree <- read.tree("data/Vascular_Plants_rooted.dated.tre")
all_plants_data <- read.csv("data/netobradley_garden_data.csv")
chemistry <- read.csv("data/neto-bradley_msc_leaf_chemistry.csv")
all_plants_aci <- read.csv("data/ecophysfit.csv")
all_plants_aci <- all_plants_aci %>% filter(Jmax != "") 
all_plants_aci$binom <- gsub('^.{16}', '', all_plants_aci$UniqueID)
all_plants_aci <- all_plants_aci %>% group_by(binom) %>% summarise(mean(Vcmax), mean(Jmax), mean(Rd))
all_plants_aci$binom <- str_to_title(all_plants_aci$binom)
all_plants_aci$binom <- gsub("_", " ", all_plants_aci$binom)

## merging data
all_plants_data$meanthick <- (all_plants_data$thick_1 + all_plants_data$thick_2 + all_plants_data$thick_3)/3
by_spp  <- group_by(all_plants_data, binomial) %>% summarise(mean_area = mean(leaf_area_imagej, na.rm = TRUE),mean_mass = mean(fresh_mass, na.rm = TRUE), mean_thickness = mean(meanthick, na.rm = TRUE))
by_spp$binomial <- gsub("_", " ", by_spp$binomial)
rownames(by_spp) <- by_spp$binomial
by_spp[130,1] <- "Tolmiea menziesii"
lookup <- lookup_table(by_spp$binomial, by_species = TRUE)
lookup$species <- rownames(lookup)
chemistry <- chemistry[1:4]
chemistry$species <- gsub("_", " ", chemistry$species)
rownames(chemistry) <- chemistry$species
by_spp <- inner_join(by_spp, chemistry,  by = c("binomial" = "species"))
by_spp <- inner_join(lookup, by_spp, by = c( "species" = "binomial"))
by_spp <- inner_join(by_spp, all_plants_aci, by = c("species" = "binom"))
by_spp <- by_spp %>% select(- sample_number)
rownames(by_spp) <- gsub(" ", "_", by_spp$species)

## make grouped data
ang <- by_spp %>% filter(group == "Angiosperms") 
rownames(ang) <- gsub(" ", "_", ang$species)
gym <- by_spp %>% filter(group == "Gymnosperms")
rownames(gym) <- gsub(" ", "_", gym$species)
pte <- by_spp %>% filter(group == "Pteridophytes")
rownames(pte) <- gsub(" ", "_", pte$species)

## cleaning up
ang_td <- treedata(all_plants_tree, ang)
gym_td <- treedata(all_plants_tree, gym)
pte_td <- treedata(all_plants_tree, pte)
##check thelypteris_novaboracensis, halocarpus bidwilii, pinus heldercichii and aristolochia tormentosa

## angiosperm data cleaning
tree <- ang_td$phy
tree <- force.ultrametric(tree)
states <- as.data.frame(ang_td$data)
traits<- states[,6:10]
traits$mean_area <- as.numeric(traits$mean_area)
traits$mean_mass <- as.numeric(traits$mean_mass)
traits$mean_thickness <- as.numeric(traits$mean_thickness)
traits$total_nitrogen_. <- as.numeric(traits$total_nitrogen_.)
traits$total_carbon_. <- as.numeric(traits$total_carbon_.)

# angiosperm ratematrix analysis
estimateTimeMCMC(data = traits, phy = tree, gen = 1000000)
#ratematrixMCMC(data = traits, phy = tree, prior = "empirical_mean", gen = 1000000, dir = "MCMC_files")
handle <- readRDS("MCMC_files/ratematrixMCMC.71327.mcmc.handle.rds")
post <- readMCMC(handle = handle, burn=0.25, thin=1000, dir = "MCMC_files")
checkConvergence(post)
logAnalyzer(handle=handle, burn=0.25, thin=1000, dir = "MCMC_files")
plotRatematrix(chain=post, colors=c("aquamarine3"))
plotRootValue(chain=post)

## gymnosperm data cleaning
tree <- gym_td$phy
tree <- force.ultrametric(tree)
states <- as.data.frame(gym_td$data)
traits<- states[,6:10]
traits$mean_area <- as.numeric(traits$mean_area)
traits$mean_mass <- as.numeric(traits$mean_mass)
traits$mean_thickness <- as.numeric(traits$mean_thickness)
traits$total_nitrogen_. <- as.numeric(traits$total_nitrogen_.)
traits$total_carbon_. <- as.numeric(traits$total_carbon_.)

# gymnosperm ratematrix analysis
estimateTimeMCMC(data = traits, phy = tree, gen = 1000000)
#ratematrixMCMC(data = traits, phy = tree, prior = "empirical_mean", gen = 1000000, dir = "MCMC_files")
handle <- readRDS("MCMC_files/ratematrixMCMC.64276.mcmc.handle.rds")
post <- readMCMC(handle = handle, burn=0.25, thin=1000, dir = "MCMC_files")
checkConvergence(post)
logAnalyzer(handle=handle, burn=0.25, thin=1000, dir = "MCMC_files")
plotRatematrix(chain=post, colors=c("aquamarine3"))
plotRootValue(chain=post)

## pteridophyte data cleaning
tree <- pte_td$phy
tree <- force.ultrametric(tree)
states <- as.data.frame(pte_td$data)
traits<- states[,6:10]
traits$mean_area <- as.numeric(traits$mean_area)
traits$mean_mass <- as.numeric(traits$mean_mass)
traits$mean_thickness <- as.numeric(traits$mean_thickness)
traits$total_nitrogen_. <- as.numeric(traits$total_nitrogen_.)
traits$total_carbon_. <- as.numeric(traits$total_carbon_.)

# pteridophyte ratematrix analysis
estimateTimeMCMC(data = traits, phy = tree, gen = 1000000)
#ratematrixMCMC(data = traits, phy = tree, prior = "empirical_mean", gen = 1000000, dir = "MCMC_files")
handle <- readRDS("MCMC_files/ratematrixMCMC.41925.mcmc.handle.rds")
post <- readMCMC(handle = handle, burn=0.25, thin=1000, dir = "MCMC_files")
checkConvergence(post)
logAnalyzer(handle=handle & handle1, burn=0.25, thin=1000, dir = "MCMC_files")
plotRatematrix(chain=post, colors=c("aquamarine3"))
plotRootValue(chain=post)




##### all the taxa in 1 plot
td <- treedata(all_plants_tree, by_spp)
tree <- td$phy
states <- as.data.frame(td$data)
trait <- states[,6:12]
trait$mean_area <- as.numeric(levels(trait$mean_area))[trait$mean_area]
trait$mean_area <- log(trait$mean_area)
trait$mean_mass <- as.numeric(levels(trait$mean_mass))[trait$mean_mass]
trait$mean_mass <- log(trait$mean_mass)
trait$mean_thickness <- as.numeric(levels(trait$mean_thickness))[trait$mean_thickness]
trait$mean_thickness <- log(trait$mean_thickness)
trait$total_nitrogen_. <- as.numeric(levels(trait$total_nitrogen_.))[trait$total_nitrogen_.]
trait$total_carbon_.<- as.numeric(levels(trait$total_carbon_.))[trait$total_carbon_.]
trait$`mean(Vcmax)`<- as.numeric(levels(trait$`mean(Vcmax)`))[trait$`mean(Vcmax)`]
trait$`mean(Jmax)`<- as.numeric(levels(trait$`mean(Jmax)`))[trait$`mean(Jmax)`]
group <- states[,4:5]
simmap <- make.simmap(tree = tree, x = group$group,model="ARD")
cols <- setNames(c("darkseagreen1","mediumseagreen","navyblue"),c("Angiosperms","Gymnosperms", "Pteridophytes"))
plotSimmap(tree=simmap, colors=cols, fsize=0.5, lwd=1)
estimateTimeMCMC(data = trait, phy = simmap, gen = 1000000)
ratematrixMCMC(data = trait, phy = simmap, prior = "empirical_mean", gen = 1000000, dir = "MCMC_files")
# mcmc chain for logged values of morpho traits
handle <- readRDS("MCMC_files/ratematrixMCMC.63338.mcmc.handle.rds")
## mcmc chain for absolute value of morpho trait
#handle <- readRDS("MCMC_files/ratematrixMCMC.69718.mcmc.handle.rds")
post <- readMCMC(handle = handle, burn=0.25, thin=1000)
checkConvergence(post)
logAnalyzer(handle=handle, burn=0.25, thin=1000)
plotRatematrix(chain=post, colors=c("#E33056","#EFC40E", "#4399DE"))
plotRootValue(chain=post)
