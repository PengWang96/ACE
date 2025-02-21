rm(list = ls())
t1 <- Sys.time()
library(GGally); library(stringr)
setwd("D:/R/factor model/real data with 6 threshold/file")
xe <- read.csv("Ebench.csv"); xv <- read.csv("Vbench.csv")
dat <- read.csv("gene_quantity_by_ACE.csv")
p <- nrow(dat)
dat <- data.frame(dat, ID = matrix(1:p,ncol = 1))

# ID_PP_not_BH <- c(2001,1432,757,776,2261,502,1947,2063,2023,444,1469,3020,1678,3251,410,1410,
#                   2926,522,3307,1343,30,562,1324,1991,1905,1957,2671,638,6,1929,1063,408,659,
#                   1009,175,2260,971,968,1209,1660,2901,2932,1456,789,2580,818,1270,3040,2457,
#                   3282,891,1231,125,3239,574,1383,1747,90,1729,1509,3319,1703,1850,1701,1532,
#                   1327,2598,2777,1601,2686,1374,1168,2701,1262,952,1977,1525,1935,1203,1988,
#                   1045,67,2116,734,783,1734,876,2263,3199,2264,1945,1039,1915,794,1075,185,
#                   368,230,1125,881,584,2058,3136,631,243,3191,2365,2222,1709,742,1344,2242,
#                   1237,2575,805,1641,3487,558)
ID_PP_not_BH <- c(2001,3064,1432,757,776,93,1383,547,2580)

# play with string, simplify the miRNA names
beginn <- str_locate(dat$miRNA, "-")
endd <- str_length(dat$miRNA)
miRNA_name <- rep('name',times = p)
for (jj in 1:p) {
  miRNA_name[jj] <- str_sub(dat$miRNA[jj], beginn[jj,1] + 1, endd[jj])
}


# data <- data.frame(t(xv - xe)); colnames(data) <- t(miRNA_name); rownames(data)<- 1:96
# data <- data[,ID_PP_not_BH]#

# data <- data.frame(t(xv)); colnames(data) <- t(miRNA_name); rownames(data)<- 1:96
# data <- data[,ID_PP_not_BH]#

data <- data.frame(t(xe)); colnames(data) <- t(miRNA_name); rownames(data)<- 1:96
data <- data[,ID_PP_not_BH]#

ggpairs(data,
        diag = "blank",
        lower = list(continuous = "smooth")
        ,title = "Correlogram for Endometrial Cancer"#Ovarian cancer    
) + 
theme(text = element_text(size = 15, hjust=0.5),
      # axis.text = element_text(angle = 45, hjust = 1, size = 10)
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size=15,hjust=0.5)
      )

# ggcorr(data, method = c("everything", "pearson"),
#        geom = "circle",
#        min_size = 1)
t2 <- Sys.time()
(t2-t1)
