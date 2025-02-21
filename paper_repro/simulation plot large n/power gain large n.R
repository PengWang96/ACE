## one sample
options(digits = 1)	
rm(list = ls())
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())

# FDR and power
dat <- read.csv("Table_one_sample.csv")

power.gain.one <- data.frame(
  Fan = (dat$power_ACE - dat$power_Fan)/dat$power_Fan * 100,
   BH = (dat$power_ACE - dat$power_BH)/dat$power_BH * 100,
  Case = dat$Case,
  var = dat$var
)


## two sample
# FDR and power
dat <- read.csv("Table_two_sample.csv")

power.gain.two <- data.frame(
  Fan = (dat$power_ACE - dat$power_Fan)/dat$power_Fan * 100,
  BH = (dat$power_ACE - dat$power_BH)/dat$power_BH * 100,
  Case = dat$Case,
  var = dat$var
)