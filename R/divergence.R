## pair picking algorithm
library(geiger)
library(dplyr)
library(ggplot2)
library(readr)

## FUNCTION FOR PLOTTING TRAIT DIV BY TAXONOMIC SPLIT
spairtraitdiv <- function(taxa, spairtime){
  library(taxonlookup)
  taxon <- lookup_table(taxa, by_species = TRUE)
  spairtime$taxgroup <- NA
  for (i in 1:nrow(spairtime)){
    tsp1 <- taxon[which(rownames(taxon) %in% spairtime$sp_1[i]),4] 
    tsp2 <- taxon[which(rownames(taxon) %in% spairtime$sp_2[i]),4] 
    taxgroup <- paste0(tsp1, tsp2)
    spairtime$taxgroup[i] <- taxgroup
  }
  spairtime$taxgroup <- gsub("GymnospermsAngiosperms", "AngiospermsGymnosperms", spairtime$taxgroup)
  spairtime$taxgroup <- gsub("PteridophytesAngiosperms", "AngiospermsPteridophytes", spairtime$taxgroup)
  spairtime$taxgroup <- gsub("PteridophytesGymnosperms", "GymnospermsPteridophytes", spairtime$taxgroup)
  ggplot(spairtime, aes(x = time_mya, y = divergence, colour = taxgroup)) + geom_point() + theme_classic()
}

## V1: USE EACH SPECIES ONCE IN A SPECIES PAIR
pick_phylo_pairs <- function(tree){
  ## label the nodes of the tree in case they are not
  tree$node.label <- c(Ntip(tree)+1:Nnode(tree))
  ## number of pairs
  n_pairs <- round(Ntip(tree)/2, digits=0)
  pairs_df <- data.frame(sp_1=rep(NA, n_pairs), sp_2=rep(NA, n_pairs))
  
  for (i in 1:(n_pairs-1)){
    print(i)
    ## pick a node at random
    tmp <- sample(tree$node.label, size=1)
    ## get descendents of nodes
    daughter <- geiger:::tips(tree, tmp)
    ## sample two of them
    spp <- sample(daughter, size=2)
    ## add them to data frame
    pairs_df[i,1] <- spp[1]
    pairs_df[i,2] <- spp[2]
    ## drop them from tree
    tree <- drop.tip(tree, tip=spp)
    ## relabel nodes
    tree$node.label <- c(Ntip(tree)+1:Nnode(tree))
  }
  ## spit out pair table
  pairs_df[n_pairs,] <- tree$tip.label
  pairs_df
}

## testing
tree <- read.tree("data/Vascular_Plants_rooted.dated.tre")
data <- read.csv("data/neto-bradley_msc_leaf_chemistry.csv")
data <- select(data, 1:4)
rownames(data) <- data$species
td <- treedata(tree, data)

my_pairs <- pick_phylo_pairs(td$phy)
datalist = list()
mrca <- mrca(phy = td$phy)
for (i in 1:length(my_pairs$sp_1)) {
  this_pair <- c(my_pairs[i,1],my_pairs[i,2])
  datalist[[i]] <- getMRCA(phy = td$phy, tip = this_pair)
}
pair_nodes <- do.call(rbind, datalist)
rownames(pair_nodes) <- pair_nodes
## extract branching times for paired taxa and add to df
ages <- as.data.frame(branching.times(td$phy))
time_mya <- ages[match(rownames(pair_nodes), rownames(ages)),]
spairtime <- cbind(time_mya, pair_nodes, my_pairs)
ggplot(spairtime, aes(time_mya)) + geom_histogram(binwidth = 50)

## for each species pair given above, calculate the variance in trait values and plot against node age of mrca
leafarea <- read.csv("data/netobradley_garden_data.csv")
leafarea <- leafarea %>% select(-c(2:7))
spairtime$divergence <- NA
for (i in 1:nrow(spairtime)) {
  sp1 <- leafarea %>% filter(binomial == spairtime$sp_1[i]) %>% summarise(mean(leaf_area_imagej))
  sp2 <- leafarea %>% filter(binomial == spairtime$sp_2[i]) %>% summarise(mean(leaf_area_imagej))
  divergence <- (log(sp1[1,]) - log(sp2[1,]))/2
  spairtime$divergence[i] <- divergence
}
spdat <- td$phy$tip.label
spairtraitdiv(taxa = spdat, spairtime = spairtime)

## V2: USE ALL POSSIBLE COMINATIONS OF SPECIES
## print all combinations of species
all_combos <- as.data.frame(expand.grid(sp_1 = spdat, sp_2 = spdat))
## remove combinations of same species
all_combos <- all_combos[which(all_combos$sp_1!=all_combos$sp_2),]


all_combos$node <- NA
ac_datalist = list()
for (i in 1:length(all_combos$sp_1)) {
  this_pair <- c(paste(all_combos[i,1]), paste(all_combos[i,2]))
  ac_datalist[[i]] <- getMRCA(phy = td$phy, tip = this_pair)
  all_combos$node[i] <- ac_datalist[[i]]
}

ac_ages <- as.data.frame(branching.times(td$phy))
colnames(ac_ages) <- "time_mya"
ac_ages$node <- rownames(ac_ages)
ac_spairtime <- merge(all_combos, ac_ages, by = c("node","node"))


ggplot(ac_spairtime, aes(time_mya)) + geom_histogram(binwidth = 50)

leafarea <- read_csv("data/netobradley_garden_data.csv")
leafarea <- leafarea %>% select(-c(2:7))
ac_spairtime$divergence <- NA
for (i in 1:nrow(ac_spairtime)) {
  sp1 <- leafarea %>% filter(binomial == ac_spairtime$sp_1[i]) %>% summarise(mean(leaf_area))
  sp2 <- leafarea %>% filter(binomial == ac_spairtime$sp_2[i]) %>% summarise(mean(leaf_area))
  ## using measure of divergence as described in Uyeda et al. 2011 
  ac_divergence <- (log(sp1[1,]) - log(sp2[1,]))/2
  ac_spairtime$divergence[i] <- as.numeric(ac_divergence)
}

ac_spairtime$unik <- paste0(ac_spairtime$time_mya, abs(ac_spairtime$divergence))
ac_spairtime <- ac_spairtime[!duplicated(ac_spairtime$unik),]
## plot by group of taxonomic divergence
spairtraitdiv(taxa = spdat, spairtime = ac_spairtime)

library(stringr)
## species pair comparison with plantecophys estimates
pep <- read_csv("data/ecophysfit.csv")
pep <- na.omit(pep)
pep$UniqueID <- str_sub(pep$UniqueID, 17)
pep$UniqueID <- str_to_title(pep$UniqueID)

ac_spairtime <- ac_spairtime[ac_spairtime$sp_1 %in% pep$UniqueID,]
ac_spairtime <- ac_spairtime[ac_spairtime$sp_2 %in% pep$UniqueID,]
for (i in 1:nrow(ac_spairtime)) {
  sp1 <- pep %>% filter(UniqueID == ac_spairtime$sp_1[i]) %>% summarise(mean(Vcmax))
  sp2 <- pep %>% filter(UniqueID == ac_spairtime$sp_2[i]) %>% summarise(mean(Vcmax))
  divergence <- (log(sp1[1,]) - log(sp2[1,]))/2
  ac_spairtime$divergence[i] <- as.numeric(divergence)
}

spdat <- td$phy$tip.label
spairtraitdiv(taxa = spdat, spairtime = ac_spairtime)

### 
## making pep2 

pep %>% filter(Vcmax > 0) %>% filter(Jmax > 0) %>% filter(UniqueID != "2019-05-29-0917_pinus_heldercichii") -> pep2

ggplot(pep2, aes(x=log(Vcmax))) + geom_histogram()
ggplot(pep2, aes(x=log(Jmax))) + geom_histogram()
library(taxonlookup)
taxa <- pep2$UniqueID
taxa[67] <- "Tolmiea_menziesii"
taxa[183] <- "Tolmiea_menziesii"
taxa
group <- vector(length=length(taxa))
ord <- vector(length=length(taxa))
for (i in 1:length(taxa)){
  g <- lookup_table(taxa[i])$group
  o <- lookup_table(taxa[i])$order
  if (length(g) == 0){
    group[i] <- NA
    ord[i] <- NA
  } else {
    group[i] <- g
    ord[i] <- o
  }
}
pep2$group <- group
pep2$order <- ord

##
## 
leafchem <- read_csv("data/neto-bradley_msc_leaf_chemistry.csv")
leafarea$leaf_area <- as.numeric(leafarea$leaf_area)
leafarea$dry_mass <- as.numeric(leafarea$dry_mass)
pep3 <- data.frame(nrow = (unique(pep2$UniqueID)))
pep3$Vcmax <- NA
pep3$Jmax <- NA
pep3$Rd <- NA
pep3$leafarea <- NA
pep3$drymass <- NA
pep3$nitr <- NA
pep3$carb <- NA
pep3$group <- NA
for (i in 1:length(leafchem$species)){
  pep3$Vcmax[i] <- pep2 %>% filter(UniqueID == leafchem$species[i]) %>% summarise(log(mean(Vcmax)))
  pep3$Jmax[i] <- pep2 %>% filter(UniqueID == leafchem$species[i]) %>% summarise(log(mean(Jmax)))
  pep3$Rd[i] <- pep2 %>% filter(UniqueID == leafchem$species[i]) %>% summarise(log(mean(Rd)))
  pep3$leafarea[i] <- leafarea %>% filter(binomial == leafchem$species[i]) %>% summarise(log(mean(leaf_area)))
  pep3$drymass[i] <- leafarea %>% filter(binomial == leafchem$species[i]) %>% summarise(log(mean(dry_mass)))
  pep3$nitr[i] <- leafchem %>% filter(species == leafchem$species[i]) %>% summarise(log(`total_nitrogen_%`))
  pep3$carb[i] <- leafchem %>% filter(species == leafchem$species[i]) %>% summarise(log(`total_carbon_%`))
  pep3$group[i] <- pep2 %>% filter(UniqueID == leafchem$species[i]) %>% summarise(first(group))
}
pep3 <- pep3[which((as.numeric(pep3$Vcmax)) != "NaN"),]
pep3$Vcmax <- as.numeric(pep3$Vcmax)
pep3$Jmax <- as.numeric(pep3$Jmax)
pep3$leafarea <- as.numeric(pep3$leafarea)
pep3$drymass <- as.numeric(pep3$drymass)
pep3$group <- as.character(pep3$group)

anova(lm(pep3$Vcmax ~ pep3$Jmax + pep3$group))
anova(lm(pep3$leafarea ~ pep3$drymass + pep3$group))
anova(lm(pep3$Vcmax ~ pep3$drymass + pep3$group))
anova(lm(pep3$Jmax ~ pep3$drymass + pep3$group))
anova(lm(pep3$leafarea ~ pep3$drymass + pep3$group))

#saveRDS(pep3, "pepthree.rds")
#gap <- lm(pep3$Vcmax ~ pep3$Jmax)
#anova(gap)

#ac_spairtime <- ac_spairtime[ac_spairtime$sp_1 %in% pep3$nrow,]
#ac_spairtime <- ac_spairtime[ac_spairtime$sp_2 %in% pep3$nrow,]
#ac_spairtime$divergence <- NA
#for (i in 1:nrow(ac_spairtime)) {
#  sp1 <- pep3 %>% filter(nrow == ac_spairtime$sp_1[i]) %>% #summarise(mean(carb))
#  sp2 <- pep3 %>% filter(nrow == ac_spairtime$sp_2[i]) %>% #summarise(mean(carb))
#  divergence <- (log(sp1[1,]) - log(sp2[1,]))/3
#  ac_spairtime$divergence[i] <- as.numeric(divergence)
#}
#ggplot(ac_spairtime, aes(x = time_mya, y = divergence)) + #geom_point()

###
library(phylolm)
td <- treedata(tree, pep3)
pep4 <- as.data.frame(td$data)
pep4$Vcmax <- as.numeric(pep4$Vcmax)
pep4$leafarea <- as.numeric(pep4$leafarea)
pep4$nitr <- as.numeric(pep4$nitr)
pep4$carb <- as.numeric(pep4$carb)
phylolm(Vcmax~leafarea+nitr+carb, data = pep4, phy = td$phy, model = "lambda")
