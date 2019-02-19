## This is code for exploring potential species sampling in a study of how morphology and metabolism are linked accross phylogeny in vascular plants
library(tidyverse)
#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)
library(plyr)

ubc <- read.csv("data/UBCBG_taxa.csv")
## 12937 rows of data
vand <- read.csv("data/VanDusen_taxa_list.csv")

head(ubc)
#taxon name // also lists genus and species seperately
head(vand)
# name (combo of g + sp + bonus:common name, usually in '...')

### things we don't want to sample:   
# known hybrids (denoted by "x " at begining, or " × " in middle for ubc, " x " in middle for vand)
# uncertain taxa (denoted by "sp.", '?' or "+")

## remove hybrids & uncertain taxa
ubc <- ubc[!grepl("x ", ubc$TaxonName),]
ubc <- ubc[!grepl("sp.", ubc$TaxonName),]
ubc <- ubc[!grepl("× ", ubc$TaxonName),]
ubc <- ubc[!grepl(" × ", ubc$TaxonName),]
ubc <- ubc[!grepl("\\?", ubc$TaxonName),]
ubc <- ubc[grepl("plant", ubc$MaterialType),]
ubc <- ubc[grepl("W", ubc$ProvenanceCode),]
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

ubc_splist <- as.vector(unique(ubc$TaxonName))
ubc_taxlist <- lookup_table(ubc_splist, by_species = TRUE)

vand_splist <- as.vector(unique(vand$NAME))
vand_taxlist <- lookup_table(vand_splist, by_species = TRUE)

## merging both lists
all_taxa_available <- rbind(vand_taxlist, ubc_taxlist)
count(all_taxa_available$family)
count(all_taxa_available$group)
count(all_taxa_available$order)
all_taxa_available$binomial <- rownames(all_taxa_available)
## 194 families represented 


## what trait data is available?
library(BIEN)
spbysp_leafn <- BIEN_trait_traitbyspecies(all_taxa_available$binomial, trait = "leaf nitrogen content per leaf area")
sub_bien <- spbysp_leafn %>% select(scrubbed_species_binomial, trait_value, project_pi)
colnames(sub_bien) <- c("taxa", "nitrogen", "source")
## 180 species
spbysp_leafa <-BIEN_trait_traitbyspecies(all_taxa_available$binomial, trait = "leaf area")
## 363 species
spbysp_leafl <-BIEN_trait_traitbyspecies(all_taxa_available$binomial, trait = "leaf life span")
## 78 species

scrubbed_tax_list <- lookup_table(unique(spbysp_leafn$scrubbed_species_binomial), by_species = TRUE)
count(scrubbed_tax_list$group) ## angio = 130 // gymno = 29 // ferns = 21
count(scrubbed_tax_list$family)

##install.packages("taxize")
#library(taxize)
#classdat <- classification(all_taxa_available$binomial, db = "tropicos", callopts = FALSE, return_id = TRUE)
print(all_taxa_available$binomial)

## figuring out where else to pull data from
try_scam <- read.csv("data/try_leafn_species_list.csv")
try_scrubbed <- try_scam[(which(try_scam$AccSpeciesName %in% all_taxa_available$binomial)),]
unique(try_scrubbed$AccSpeciesName)
summary(try_scrubbed$Dataset)

## largest datasets
## reich-oleksyn global leaf n, p database ## downloaded 2019Feb04
## topic ## not open access 2019Feb05
## plant volatiles ## n data is not open access 2019Feb05
## global 15N database ## not open access 2019Feb05
## catalonian mediterranean forest trait db ## downloaded 2019Feb05

reich_dat <- read.csv("data/reich_oleksyn_dataset.csv")
reich_dat <- reich_dat[which(reich_dat$Species %in% all_taxa_available$binomial),]
sub_reich <- reich_dat %>% select(Species, N..mg.g, References)
colnames(sub_reich) <- c("taxa", "nitrogen", "source")
## 153 taxa in reich dataset for N content

brot <- read.csv("data/BROT2_dat.csv")
brot <- brot %>% filter(Trait == "LNCm")
brot <- brot[(which(brot$Taxon %in% all_taxa_available$binomial)),]
sub_brot <- brot %>% select(Taxon, Data, SourceID)
colnames(sub_brot) <- c("taxa", "nitrogen", "source")

glopnet <- read.csv("data/glopnet.csv")
glopnet <- glopnet[which(glopnet$Species %in% all_taxa_available$binomial),]
sub_glop <- glopnet %>% select(Species, log.Nmass, Dataset)
colnames(sub_glop) <- c("taxa", "nitrogen", "source")

all_data_available <- rbind(sub_glop, sub_reich, sub_brot, sub_bien)
## data available for 272 taxa
## nitrogen values not in same units - check for log transformation in glopnet dataset

taxonlist <- as.character(unique(all_data_available$taxa))
taxonlist <- lookup_table(taxonlist, by_species = TRUE)

## visualizing species sample spread 
library(ggplot2)
ggplot(data = taxonlist, aes(x = group, fill = order)) + geom_histogram(stat = "count") + theme_minimal() + xlab("Group") + ylab("Number of Species")
ggplot(data = taxonlist, aes(x = family, fill = family)) + geom_histogram(stat = "count") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),) 


## getting other traits
all_data_available$taxa <- as.character(all_data_available$taxa)
traits <- BIEN_trait_list()
traits <- traits$trait_name
traits <- as.character(traits)

ll <- BIEN_trait_traitbyspecies(all_data_available$taxa, trait = traits)
ll_a <- ll %>% filter(ll$trait_name == "leaf life span")
ll_a$scrubbed_species_binomial <- as_factor(ll_a$scrubbed_species_binomial)
## trait value needs to be numeric
ll_a$trait_value <- as.numeric(ll_a$trait_value)

