#optimization data analysis 
#KG

setwd("/Users/Katie/Desktop/R files from work/270717 Optimisation/")
getwd()

require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)

load("ATT93447.RData")

#For each antigen, use the factor level


