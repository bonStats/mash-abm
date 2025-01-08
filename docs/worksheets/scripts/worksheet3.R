
library(BMS)

library(readr)
library(dplyr)
library(ggplot2)


josh_data_dir <- function(fl){
  paste0("/Users/jbon/github/mash-abm/data/", fl)
}

deathrate <- read_csv(file = josh_data_dir("deathrate2.csv"), col_types = cols(col_skip(), col_guess()))


# Bayesian model averaging with BMS

reg <- bms(data.frame(deathrate), burn=1e4, iter=5e4)
coef(reg)
image(reg)
topmodels.bma(reg)[, 1:5]

