library(dplyr)

correct_racir <- function(empty, file){

  raci <- read_li68(empty)

# make timescale numeric
raci$elapsed <- as.numeric(raci$elapsed)

# remove first and last minute of logged data (non linear section)
raci <- raci %>% 
  arrange(elapsed) %>% 
  filter(elapsed > 60) %>%
  filter(elapsed < 600)

## visual prompt to check that only the linear parts of the ramp are kept
plot(raci$CO2_r, raci$A)

## estimate parameters for correction of curve data
model <- lm(raci$A ~ poly(raci$CO2_r,2, raw = TRUE))
b_cor <- as.numeric(coefficients(model)[2])
a_cor <- as.numeric(coefficients(model)[3])
c_cor <- as.numeric(coefficients(model)[1])


## read in the real data          
d_aci <- read_li68(file)

# make timescale numeric
d_aci$elapsed <- as.numeric(d_aci$elapsed)

# trim first and last minute of logged data (non linear section of ramp)
d_aci <- d_aci %>% 
  arrange(elapsed) %>% 
  filter(elapsed > 60) %>%
  filter(elapsed < 600)

## visual prompt to check that only linear parts of the ramp are kept
plot(d_aci$CO2_r, d_aci$A)    

## takes parameters from fitted polynomial and calculate correction value
d_aci$correction <- ((a_cor*(d_aci$CO2_r^2))+(b_cor*d_aci$CO2_r)+c_cor)

## adds a corrected A value
d_aci$Aleaf <- d_aci$A-d_aci$correction


## need to add a corrected Ci value!!

return(d_aci)

}

##testing
therealdata <- correct_racir(empty = "data/practice_data/2019-04-10-1050_racir_empty", file = "data/practice_data/2019-04-10-1109_racir_prayerplant")


## should add this in to the above code
## therealdata$Aleaf <- therealdata$A-therealdata$correction
