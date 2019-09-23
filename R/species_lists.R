## Code for determining species sampling strategy 

library(tidyverse)
#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)
library(plyr)

# full garden lists
ubc <- read.csv("data/species_list/UBCBG_taxa.csv")
vand <- read.csv("data/species_list/VanDusen_taxa_list.csv")

### things we don't want to sample:   
# known hybrids (often denoted by "x ")
# uncertain taxa (denoted by "sp.", '?' or "+")

## remove hybrids & uncertain taxa
ubc <- ubc[!grepl("x ", ubc$TaxonName),]
ubc <- ubc[!grepl("sp.", ubc$TaxonName),]
ubc <- ubc[!grepl("× ", ubc$TaxonName),]
ubc <- ubc[!grepl(" × ", ubc$TaxonName),]
ubc <- ubc[!grepl("\\?", ubc$TaxonName),]
ubc <- ubc[grepl("plant", ubc$MaterialType),]
ubc <- ubc[!grepl("var", ubc$TaxonName),]
ubc <- ubc[!grepl("\\'", ubc$TaxonName),]

vand <- vand[!grepl("x ", vand$NAME),]
vand <- vand[!grepl("sp.", vand$NAME),]
vand <- vand[!grepl("× ", vand$NAME),]
vand <- vand[!grepl(" × ", vand$NAME),]
vand <- vand[!grepl("\\+", vand$NAME),]
vand <- vand[!grepl("var.", vand$NAME),]
vand <- vand[!grepl("f.", vand$NAME),]
vand <- vand[!grepl("Group", vand$NAME),]
vand <- vand[!grepl("\\'", vand$NAME),]
vand <- vand[!grepl("\\(", vand$NAME),]
vand <- vand[!grepl("garden origin", vand$NATIVITY),]
vand <- vand[!grepl("Garden Origin", vand$NATIVITY),]

## get full taxonomic listing of the remaining species
splist <- c(as.vector(unique(ubc$TaxonName)),as.vector(unique(vand$NAME)))
splist <- unique(splist)
splist <- lookup_table(splist, by_species = TRUE)

## summary stats
splist$binomial <- rownames(splist)
count(splist$family)
count(splist$group)
count(splist$order)

## check for species overlap with Zanne et al. tree 
library(ape)
library(geiger)
phy <- read.tree(file = "data/Vascular_Plants_rooted.dated.tre")
dat <- as.matrix(splist)
splist$binomial <- gsub(" ", "_", splist$binomial)
rownames(dat) <- splist$binomial
td <- treedata(phy, dat)
all_taxa_available_phy <- rownames(td$data)

all_taxa_available_phy <- lookup_table(all_taxa_available_phy, by_species = TRUE)
all_taxa_available_phy$binomial <- rownames(all_taxa_available_phy)

# add source garden id
ubc$TaxonName <- gsub(" ", "_", ubc$TaxonName)
for (i in 1:nrow(all_taxa_available_phy)) { 
  if(all_taxa_available_phy$binomial[i] %in% ubc$TaxonName) {
    all_taxa_available_phy$id[i] <- "ubc"
  } 
  else { 
    all_taxa_available_phy$id[i] <- "vand"
  }
}

## summary stats
count(all_taxa_available_phy$group) ## angio = 1594 // gymno = 110 // ferns = 38
count(all_taxa_available_phy$family)

## generating sample
# garden_sample <- data_frame()
# ufam <- unique(all_taxa_available_phy$family)
# gsample <- data.frame(matrix(nrow = 206, ncol = 6))
# colnames(gsample) <- colnames(all_taxa_available_phy)
#for (i in 1:length(ufam)) {
#  garden_sample <- all_taxa_available_phy %>% 
#  dplyr::filter(family == ufam[i]) 
#  gsample[i, ] <- (sample_n(garden_sample,size = 1))
# }

#gsample[180:188,] <- sample_n((all_taxa_available_phy %>% 
  #          dplyr::filter(family == "Ericaceae")) , size = 9)
#gsample[189:197,] <- sample_n((all_taxa_available_phy %>% 
  #          dplyr::filter(family == "Pinaceae")) , size = 9)
#gsample[198:206,] <- sample_n((all_taxa_available_phy %>% 
 #           dplyr::filter(family == "Dryopteridaceae")) , size = 9)

#gsample$binomial <- gsub(x = gsample$binomial, "_", " ")
#vdsample <- gsample %>% dplyr::filter(id == "vand")
#ubcsample <- gsample %>% dplyr::filter(id == "ubc")
#write.csv(ubcsample, file = "netobradley_sample_ubc.csv", row.names = FALSE)
#write.csv(vdsample, file = "netobradley_sample_vd.csv", row.names = FALSE)
#write.csv(gsample, file = "netobradley_sample.csv")


all_taxa_available_phy %>% 
  dplyr::filter(family == "Montiaceae") %>%
  sample_n(size = 1)


### replacing sample

garden_sample <- data_frame()
ufam <- replace$family
gsample <- data.frame(matrix(nrow = 29, ncol = 6))
all_ubc <- all_taxa_available_phy[all_taxa_available_phy$id == "ubc",]
all_ubc <- as.data.frame(all_ubc)
colnames(gsample) <- colnames(all_ubc)
for (i in 1:29) {
 garden_sample <- all_ubc %>% dplyr::filter(all_ubc$family == 'Ericaceae')
sample_n(garden_sample,size = 1, replace = FALSE)


garden_sample <- data_frame()
ufam <- replace$family
gsample <- data.frame(matrix(nrow = 29, ncol = 6))
all_vd <- all_taxa_available_phy[all_taxa_available_phy$id == "vand",]
all_vd <- as.data.frame(all_vd)
colnames(gsample) <- colnames(all_vd)
for (i in 1:29) {
  garden_sample <- all_ubc %>% dplyr::filter(all_ubc$family == 'Ericaceae')
  sample_n(garden_sample,size = 1, replace = FALSE)
  