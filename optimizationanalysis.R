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

# set negative normalized values to zero

rep1.matrix <- as.matrix(trans.norm.rep1)

for (i in 1:length(rep1.matrix))
{
  if (is.na(rep1.matrix[[i]]) | rep1.matrix[[i]] > 0) 
  {
    i = i+1
  } else if(rep1.matrix[[i]] < 0) 
  {
    rep1.matrix[[i]] <- 0
  }
}

rep2.matrix <- as.matrix(trans.norm.rep2)

for (i in 1:length(rep2.matrix))
{
  if (is.na(rep2.matrix[[i]]) | rep2.matrix[[i]] > 0) 
  {
    i = i+1
  } else if(rep2.matrix[[i]] < 0) 
  {
    rep2.matrix[[i]] <- 0
  }
}


#average replicates - replicates are stored in trans.norm.rep1 and trans.norm.rep2
### Average duplicates, if the data has technical replicates in the form of 2 blocks / subarray
  
  trans.norm.avg <- matrix(nrow = nrow(trans.norm.rep1), ncol = ncol(trans.norm.rep1))
  colnames(trans.norm.avg) = colnames(trans.norm.rep1)
  rownames(trans.norm.avg) = rownames(trans.norm.rep1)
  
  trans.norm.avg <- log2((2^rep1.matrix + 2^rep2.matrix)/2)
  

### Check for deviant technical replicates, automatically exclude (set to NA)
# Use Patrick's formula for ELISA to compare replicates within one array
# if rep1 or rep2 is more than 1.5 times rep2 or rep1, respectively, exclude that pair
# Also, can redo this using the subsetted matrices (rep1 and rep2) and it should be shorter

 for(k in 1:ncol(rep1.matrix)){
   for(j in 1:nrow(rep1.matrix)){
     if(is.na(rep1.matrix[j,k]) | is.na(rep2.matrix[j,k]) | rep1.matrix[j,k]<0.2 | rep2.matrix[j,k]<0.2){
       j+1
     } else if (rep1.matrix[j,k] > (log2(1.5) + rep2.matrix[j,k]) | (rep2.matrix[j,k] > (log2(1.5) + rep1.matrix[j,k])) == TRUE) 
     {
       trans.norm.avg[j,k] <- NA
     }
   }
 }
remove(j,k)
  
#subtract GST - subtract GST at the same dilution as the other antigen was diluted
#need to deal with NAs this time

# ratio of positive to negative - subtract the negative log2 value for each condition



