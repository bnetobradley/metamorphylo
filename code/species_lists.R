## This is code for exploring potential species sampling in a study of how morphology and metabolism are linked accross phylogeny in vascular plants
library(tidyverse)
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)

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
## 194 families represented 
