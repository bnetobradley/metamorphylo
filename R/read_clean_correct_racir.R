## sourcing from read licor 6800 script for appropriate file formatting
source('R/read_licor.R')
source('R/recalculate_licor_racir.R')
source('R/correct_racir.R')
library(ggplot2)
library(dplyr)

## read all txt files and merge into one dataframe: original file source denoted by id column
nm <- list.files(path="data/licor_data/txt", full.names = TRUE)
lfread <- do.call(rbind, lapply(nm, read_li68))

## new column in licor dataframe with date in yyyy-mm-dd format 
# add column with true date
lfread$id <- gsub("data/licor_data/txt/","", lfread$id)
lfread$date_collected <- substr(lfread$id, 0, 10)
lfread <- lfread %>% select(-X147)

## new column in licor dataframe with species name in genus_species format
splicor <- c()
for (i in 1:length(lfread$id)) {
  splicor[i] <- (word(lfread$id[i], 1:2, -1, sep = fixed("_"))[2])
}
lfread$binomial <- splicor

##ADD LEAF AREA IN CHAMBER VALUE TO RACIR DATASET
## read in leaf area in chamber data 
leafarea <- read.csv("data/netobradley_garden_data.csv")

## new column in leaf area dataframe with date in yyyy-mm-dd format 
format(Sys.Date(), "%Y%b%d")
leafarea$date_collected <- as.Date(leafarea$date_collected, "%Y%b%d")

## leaf area binomial in lower case in genus_species format
leafarea$binomial <- tolower(leafarea$binomial)

## create a 'unique' matching column with binomial and date
leafarea$uniquespdate <- paste(leafarea$binomial, leafarea$date_collected)
lfread$uniquespdate <- paste(lfread$binomial, lfread$date_collected)
## for each leaf sampled: add the correct leaf area in chamber to the licor file
## if unique codes match, paste chamber leaf area in col

## confirm species names spelling matches in both files
unique(lfread$binomial[((lfread$binomial %in% leafarea$binomial) == FALSE)])

## add column with mean projected leaf area in chamber to licor files
leafarea$area_correction <- ifelse(leafarea$chamber_area=="full", print(leafarea$chamber_size_cm2), print((leafarea$circle_1 + leafarea$circle_2 + leafarea$circle_3)/3))

# merge leaf area correction with licor files
check <- leafarea %>% select(area_correction,uniquespdate)
lfread <- full_join(lfread, check, by = "uniquespdate" )

# replace S (leaf area value) if not empty
 for( i in 1:nrow(lfread)){
   ifelse(!is.na(lfread$area_correction[i]) == TRUE,
          lfread$S[i] <- lfread$area_correction[i],
          lfread$S[i] <- lfread$S[i]  )
 } 

lfread <- lfread %>% filter(S < 7)
newlfread <- recalc_licor68(lfread)
        
## make a list of unique measurements and the closest empty curve (in time)
## unique curve lists and new df for data and empty curves
empty_curves <- newlfread[which(newlfread$binomial == "empty_racir"),]
data_curves <- newlfread[which(newlfread$binomial != "empty_racir"),]
unq <- unique(data_curves$id)
emp_unq <- unique(empty_curves$id)
outs <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(outs) <-  c("numb", "unique_id")

## match each data curve to empty curve by proximity in time
## outs - list of paired curves
for (i in 1:length(unq)) { 
  temp_curves <- data_curves %>% filter(id == unq[i])
  outp <- first(which(abs(empty_curves$TIME - temp_curves$TIME[1]) == min(abs(empty_curves$TIME - temp_curves$TIME[1]), na.rm = TRUE)))
  out <- empty_curves[outp, 'id']
  ucode <- print(unq[i]) 
  outs[i,1] <- out
  outs[i,2] <- ucode
}

# save cleaned data files and empty/data matching list
#write.csv(outs, "emptydatamatch.csv")
#write.csv(newlfread, "areacorrectracirfiles.csv")