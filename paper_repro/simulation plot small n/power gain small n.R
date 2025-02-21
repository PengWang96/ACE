## one sample
setwd("D:/R/factor model20230625/simulation plot small n")
options(digits = 1)	
rm(list = ls())
# FDR and power
dat <- read.csv("Table_one_sample.csv")

power.gain.one <- data.frame(
  Fan = (dat$power_ACE - dat$power_Fan)/dat$power_Fan * 100,
  BH = (dat$power_ACE - dat$power_BH)/dat$power_BH * 100,
  Case = dat$Case,
  var = dat$var
)


## two sample
options(digits = 1)	
# FDR and power
dat <- read.csv("Table_two_sample.csv")

power.gain.two <- data.frame(
  Fan = (dat$power_ACE - dat$power_Fan)/dat$power_Fan * 100,
  BH = (dat$power_ACE - dat$power_BH)/dat$power_BH * 100,
  Case = dat$Case,
  var = dat$var
)