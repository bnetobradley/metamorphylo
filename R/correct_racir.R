library(dplyr)

source(file = "R/read_licor.R")

correct_racir <- function(empty, file) {
  raci <- empty

# remove first and last minute of logged data (non linear section)
#raci <- raci %>% 
 # arrange(elapsed) %>% 
  #filter(elapsed > 60) %>%
  #filter(elapsed < 600)

## visual prompt to check that only the linear parts of the ramp are kept
plot(raci$CO2_r, raci$A, xlim = c(0,1100), ylim = c(-100,50))

## read in the real data          
d_aci <- file

# trim first and last minute of logged data (non linear section of ramp)
#d_aci <- d_aci %>% 
 # arrange(elapsed) %>% 
  #filter(elapsed > 60) %>%
  #filter(elapsed < 600)

## visual prompt to check that only linear parts of the ramp are kept
#plot(d_aci$CO2_r, d_aci$A, xlim = c(0,1100), ylim = c(-100,50))    
  
## estimate parameters for correction of curve data
cal1st <- lm(A ~ CO2_r, data = raci)
cal2nd <- lm(A ~ poly(CO2_r, 2), data = raci)
cal3rd <- lm(A ~ poly(CO2_r, 3), data = raci)
cal4th <- lm(A ~ poly(CO2_r, 4), data = raci)
cal5th <- lm(A ~ poly(CO2_r, 5), data = raci)
bics <- BIC(cal1st, cal2nd, cal3rd, cal4th, cal5th)
best <- noquote(rownames(bics)[bics$BIC == min(bics$BIC)])

## takes parameters from fitted polynomial and

ifelse(best == "cal5th",
  d_aci$Acor <- (d_aci$A - predict(cal5th)),
ifelse(best == "cal4th",
  d_aci$Acor <- (d_aci$A - predict(cal4th)),
ifelse(best == "cal3rd",
  d_aci$Acor <- (d_aci$A - predict(cal3rd)),
ifelse(best == "cal2nd",
  d_aci$Acor <- (d_aci$A - predict(cal2nd)),
ifelse(best == "cal1st",
  d_aci$Acor <- (d_aci$A - predict(cal1st)))))))

d_aci$Cicor <- (((d_aci$gtc - d_aci$E/2)* d_aci$Ca - d_aci$Acor)/(d_aci$gtc + d_aci$E/2))

plot(d_aci$Cicor, d_aci$Acor, xlim = c(0,1100), ylim = c(-10,50))

return(d_aci)

}

##testing
#therealdata <- correct_racir(empty = raccal, file = racdata)
