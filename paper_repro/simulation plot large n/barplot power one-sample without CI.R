rm(list = ls())
library(scales)
library(tidyverse)
library(viridis)
library(gcookbook)
library(gridExtra)
library(dplyr)
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
dat1 <- dat[which(dat$var == "Homoscedasticity"), ]
dat2 <- dat[which(dat$var == "Heteroscedasticity"), ]
##############  Homoscedasticity
# FDR under Homoscedasticity
data1.FDR <- data.frame(
  individual = dat1$Case,
  group = as.factor(rep(c('C1', 'C2', 'C3','C4'), times = nrow(dat1)/4)),
  ACE = dat1$FDR_ACE,
  PFA = dat1$FDR_Fan,
   BH = dat1$FDR_BH,
  `e-BH` = dat1$FDR_e_BH,
  BY = dat1$FDR_BY
)
data1.FDR <- data1.FDR %>% gather(key = "Method", value="value", -c(1,2))
data1.FDR$Test <- "Homoscedasticity"

# power under Homoscedasticity
data1.power <- data.frame(
  individual = dat1$Case,
  group = as.factor(rep(c('C1', 'C2', 'C3','C4'), times = nrow(dat1)/4)),
  ACE = dat1$power_ACE,
  PFA = dat1$power_Fan,
  BH = dat1$power_BH,
  `e-BH` = dat1$power_e_BH,
  BY = dat1$power_BY
)
data1.power <- data1.power %>% gather(key = "Method", value="value", -c(1,2))
data1.power$Test <- "Homoscedasticity"

############### Heteroscedasticity
# FDR under Heteroscedasticity
data2.FDR <- data.frame(
  individual = dat2$Case,
  group = as.factor(rep(c('C1', 'C2', 'C3','C4'), times = nrow(dat2)/4)),
  ACE = dat2$FDR_ACE,
  PFA = dat2$FDR_Fan,
  BH = dat2$FDR_BH,
  `e-BH` = dat2$FDR_e_BH,
  BY = dat2$FDR_BY
)
data2.FDR <- data2.FDR %>% gather(key = "Method", value="value", -c(1,2))
data2.FDR$Test <- "Heteroscedasticity"

# power under Heteroscedasticity
data2.power <- data.frame(
  individual = dat2$Case,
  group = as.factor(rep(c('C1', 'C2', 'C3','C4'), times = nrow(dat2)/4)),
  ACE = dat2$power_ACE,
  PFA = dat2$power_Fan,
  BH = dat2$power_BH,
  `e-BH` = dat2$power_e_BH,
  BY = dat2$power_BY
)
data2.power <- data2.power %>% gather(key = "Method", value="value", -c(1,2))
data2.power$Test <- "Heteroscedasticity"

data1.FDR$Method <- gsub("e.BH", "e-BH", data1.FDR$Method)
data2.FDR$Method <- gsub("e.BH", "e-BH", data2.FDR$Method)
data1.power$Method <- gsub("e.BH", "e-BH", data1.power$Method)
data2.power$Method <- gsub("e.BH", "e-BH", data2.power$Method)

# add expression
f_names <- list('C1' = expression(paste(pi[1] == 0.2, "  ", B != 0)), 
                'C2' = expression(paste(pi[1] == 0.2,"  ", B == 0)), 
                'C3' = expression(paste(pi[1] == 0.05,"  ", B != 0)), 
                'C4' = expression(paste(pi[1] == 0.05,"  ", B == 0)))
f_labeller <- function(variable, value){return(f_names[value])}

# Make the plot
p1.FDR <- ggplot(data1.FDR, aes(x=as.factor(individual), y=value, fill=Method)) +      
  geom_bar(stat="identity", width=0.7, position = position_dodge(), color="black") +
  geom_hline(aes(yintercept=0.05), linetype = "dashed", size = 1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.06)) +
  scale_fill_manual(values = c("#F94141", "#1868B2", "#018A67", "#F98F34", "#6A4C93")) +
  theme_bw(base_size = 14) + 
  facet_wrap( . ~ group, ncol = 4, labeller = f_labeller) +
  theme(legend.position = "none", 
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title=element_blank()) +
  labs(title = "Homoscedasticity", x = NULL, y = "FDR")
p1.FDR


p2.FDR <- ggplot(data2.FDR, aes(x=as.factor(individual), y=value, fill=Method)) +      
  geom_bar(stat="identity", width=0.7, position="dodge", color="black") +
  geom_hline(aes(yintercept=0.05), linetype = "dashed", size = 1) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 0.06)) +
  scale_fill_manual(values = c("#F94141", "#1868B2", "#018A67", "#F98F34", "#6A4C93")) +
  theme_bw(base_size = 15) + 
  facet_wrap(. ~ group, ncol = 4, labeller = f_labeller) +
  theme(legend.position = "none", 
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title=element_blank()) +
  labs(title = "Heteroscedasticity", x = "(a) FDR control", y = "FDR")
p2.FDR

p1.power <- ggplot(data1.power, aes(x=as.factor(individual), y=value, fill=Method)) +      
  geom_bar(stat="identity",width=0.7, position="dodge", color="black") +
  scale_fill_manual(values = c("#F94141", "#1868B2", "#018A67", "#F98F34", "#6A4C93")) +
  theme_bw(base_size = 15) + 
  facet_wrap(. ~ group, ncol = 4, labeller = f_labeller) +
  theme(legend.position = "none", 
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title=element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0.25, 0.97), oob = rescale_none) + 
  labs(title = "Homoscedasticity", x = NULL, y = "Power")
p1.power

p2.power <- ggplot(data2.power, aes(x=as.factor(individual), y=value, fill=Method)) +      
  geom_bar(stat="identity",width=0.7, position="dodge", color="black") +
  scale_fill_manual(values = c("#F94141", "#1868B2", "#018A67", "#F98F34", "#6A4C93")) +
  labs(title = "Heteroscedasticity", x = "(b) Power comparison", y = "Power") + 
  theme_bw(base_size = 15) + 
  facet_wrap(. ~ group, ncol = 4, labeller = f_labeller) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 12),
        legend.key.width = unit(55, "pt"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title=element_blank(),
        legend.text = element_text(face="italic", size = 15)
        ) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.82), oob = rescale_none)
p2.power

p <- grid.arrange(arrangeGrob(p1.FDR, p2.FDR, p1.power, p2.power, ncol=1,
                         heights=c(1, 1.1, 1.2, 1.7)))

