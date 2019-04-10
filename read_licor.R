## Adapted from code written in 2018 by Richard Telford and Barbara Neto-Bradley for Temperature response data at PFTC-4
## adapted for RACiR data on LI-6800 by BNB in 2019

library(tidyverse)
read_li68 <- function(file, debug = FALSE){
  
  if(debug){
    print(file)
  }
  ## reading in raw licor data file
  f <- read_lines(file = file) 
  #find data_row
  data_row <- grep("\\[Data]", f)
  #remove everything before data rows (also remove 1 row below data header)
  f <- f[-(1:(data_row + 1))]  
  f <- f[-(2)]
  

  #remove everything below end header row 
  head_row <- grep("\\[Header]", f)
  if (length(head_row) > 0){
  f <- f[-((head_row):(length(f)))]  
  }
  
  dat <- read_delim(paste(f, collapse = "\n"), delim = "\t")
  
  return(dat)
}
