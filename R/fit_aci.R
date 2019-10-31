#install.packages("plantecophys")
library(plantecophys)
library(bayCi)
library(units)
library(magrittr)
source(file = "R/correct_racir.R")

## Fit an A-Ci curve with plantecophy package
## racir data needs to be pre-corrected for this

# read in area corrected data file
areacorracir <- read_csv("~/Desktop/areacorrectracirfiles.csv")

## read in dataframe with data curves and empty curves matched by time
empdatmatch <- read.csv("data/emptydatamatch.csv")
empdatmatch <- empdatmatch %>% select(numb, unique_id)

## get rid of empty and data curves that have less than 330 observations 
## start with 170349 observations
unq <- unique(areacorracir$id)
removefile <- c()
for (i in 1:length(unq)) {
  data <- areacorracir %>% filter(id == unq[i])
  ifelse((nrow(data) >= 330),
         removefile <- NA, 
         removefile <- print(unq[i]))
  ifelse((is.na(removefile)),
         areacorracir <- areacorracir , 
    areacorracir <- areacorracir[-which(areacorracir$id == unq[i]),])
}
## end up with 168340 observations


## Fit A-Ci curves using ecophys package
pep_fit <- data.frame(matrix(nrow = 430, ncol = 4))
colnames(pep_fit) <- c("Vcmax", "Jmax", "Rd", "UniqueID")
for (i in 1:nrow(empdatmatch)){
  tryCatch({
    empty <- areacorracir %>% 
      filter(id == empdatmatch$numb[i]) %>% 
      mutate(
        A  = set_units( A,  umol / m^2 / s ),
        Cr = set_units( CO2_r, umol / mol     ),
        time = set_units(time, s)
      ) %>% 
      bayCi:::prepare_empty() %>%
      dplyr::mutate_if(~inherits(.x, "units"), units::drop_units) 
    empty %<>% filter(empty$use == TRUE)
    notempty <- areacorracir %>% filter(id == empdatmatch$unique_id[i]) %>% filter(obs %in% empty$obs)
    corrected <- correct_racir(empty = empty, file = notempty)
    fitecophys <- fitaci(data = corrected, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", PPFD = "Q"))
    # open file for saving plot
    mypath <- paste("data/ecophysplots/", empdatmatch$unique_id[i], ".jpg", sep = "")
    jpeg(file=mypath)
    plot(fitecophys, sub = empdatmatch$unique_id[i])
    dev.off()
    pep_fit[i, 1] <- (fitecophys$pars)[1,1]
    pep_fit[i, 2] <- (fitecophys$pars)[2,1]
    pep_fit[i, 3] <- (fitecophys$pars)[3,1]
    pep_fit[i, 4] <- paste(empdatmatch$unique_id[i])
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#write_csv(pep_fit, "data/ecophysfit.csv")


######
######


## Plot data for which ecophys package won't fit curve
for (i in 1:nrow(not_fitting)) {
  # empty_nf <- areacorracir %>% filter(id == not_fitting$empty[i])
  #  data_nf <- areacorracir %>% filter(id == not_fitting$data[i])
  # ggplot(data = data_nf, aes(x=CO2_r, y=A)) + geom_point() + geom_point(data = empty_nf, aes(x=CO2_r, y=A), colour = "red") + ylim(0,20)
  
  empty_nf <- areacorracir %>% 
    filter(id == not_fitting$empty[i]) %>% 
    mutate(
      A  = set_units( A,  umol / m^2 / s ),
      Cr = set_units( CO2_r, umol / mol     ),
      time = set_units(time, s)
    ) %>% 
    bayCi:::prepare_empty() %>%
    dplyr::mutate_if(~inherits(.x, "units"), units::drop_units) 
  empty_nf %<>% filter(empty_nf$use == TRUE)
  notempty_nf <- areacorracir %>% filter(id == not_fitting$data[i]) %>% filter(obs %in% empty_nf$obs)
  corrected <- correct_racir(empty = empty_nf, file = notempty_nf)
  mypath <- paste("data/ecophysplots/failing/", not_fitting$data[i], ".jpg", sep = "")
  jpeg(file=mypath)
  print(ggplot(data = corrected, aes(x=Cicor, y=Acor)) + geom_point(colour = "purple") + labs(subtitle = paste(not_fitting$data[i], not_fitting$empty[i], sep = " & ") ) + geom_point(data = empty_nf, aes(x = CO2_r, y = A), colour = "red") + geom_point(data = notempty_nf, aes(x = CO2_r, y = A), colour = "blue"))
  dev.off()
  
}

areacorracir %>% filter(id == not_fitting[i])