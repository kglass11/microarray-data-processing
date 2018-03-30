#optimization data analysis 
#KG

setwd("/Users/Katie/Desktop/R files from work/270717 Optimisation/")
getwd()

# install.packages("lme4")
library(lme4)

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

#import target metadata file
target_file <- "Opto Target Metadata.csv" 
target_meta.df <- read.csv(target_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

# set negative normalized values to zero - replicates are stored in trans.norm.rep1 and trans.norm.rep2
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

### Average duplicates
  trans.norm.avg <- matrix(nrow = nrow(trans.norm.rep1), ncol = ncol(trans.norm.rep1))
  colnames(trans.norm.avg) = colnames(trans.norm.rep1)
  rownames(trans.norm.avg) = rownames(trans.norm.rep1)
  
  trans.norm.avg <- log2((2^rep1.matrix + 2^rep2.matrix)/2)

### Check for deviant technical replicates, automatically exclude (set to NA)
# Use Patrick's formula for ELISA to compare replicates within one array
# if rep1 or rep2 is more than 1.5 times rep2 or rep1, respectively, exclude that pair
# Added: if either rep1 or rep2 is less than 0.2, then don't apply the rule
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
  # subtract GST from all except AMA1; that means MSP1-19, Hyp2, GEXP18,

#merge target metadata with trans.norm.avg (transposed)
target.df <- merge(target_meta.df, t(trans.norm.avg), by.x = "Name", by.y = "row.names", sort = FALSE)

#Extract GST elements for all print buffers and concentrations
GST <- c(grep("GST", target.df$Name))
GST.df <- target.df[GST,]

#function for subtracting GST from a concentration and print buffer
subtractGST <- function(concentration, print_buffer){
  #isolate data for that concentration and print buffer
  con.df <- filter(target.df, Concentration == concentration, Print_Buffer == print_buffer, Expression_Tag == "GST")
  con.df <- tibble::column_to_rownames(con.df, var="Name")
  con.df <- con.df[,(ncol(con.df) - ncol(t(trans.norm.avg))+1):ncol(con.df)]
  
  GSTsub.df <- filter(GST.df, Concentration == concentration, Print_Buffer == print_buffer)
  GSTsub.df <- tibble::column_to_rownames(GSTsub.df, var="Name")
  GSTsub.df <- GSTsub.df[,(ncol(GSTsub.df) - ncol(t(trans.norm.avg))+1):ncol(GSTsub.df)]
  
  con2.df <- data.frame(matrix(0, nrow = nrow(con.df), ncol = ncol(con.df)))
  rownames(con2.df) <- rownames(con.df)
  colnames(con2.df) <- colnames(con.df)
  
  #subtract GST!
  for(b in 1:ncol(con2.df)){
    for(a in 1:nrow(con2.df)){
      #if GST value is NA or main value is NA, then cannot do subtraction and cannot compare with other data
      #so set those values to NA
      if(is.na(GSTsub.df[1,b])|is.na(con.df[a,b])){
        con2.df[a,b] <- NA
        # when the GST value is positive only, subtract GST value, 
        #otherwise want to leave as whatever the value was before (because GST was at or below buffer mean for that sample)
      } else if (GSTsub.df[1,b] > 0){
        #calculate difference in original MFI form (not log2)
        con2.df[a,b] <- 2^con.df[a,b] - 2^GSTsub.df[1,b]
        #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
        if (con2.df[a,b] > 0) {
          con2.df[a,b] <- log2(con2.df[a,b])
          #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
          #so also need to set those negatives to 0 again.
          if(con2.df[a,b] < 0){
            con2.df[a,b] <- 0
          }
        } else { 
          con2.df[a,b] <- 0
        }
      } else {
        con2.df[a,b] <- con.df[a,b]
      }
    }
  }
  remove(a,b)
  
  #return GST-subtracted data
  return(con2.df)
  
}

#Subtract GST from each print buffer and concentration
#I realize it would be easier to use a loop or something here but I was tired
GST100B1.df <- subtractGST(100,1)
GST100B2.df <- subtractGST(100,2)
GST100B3.df <- subtractGST(100,3)

GST25B1.df <- subtractGST(25,1)
GST25B2.df <- subtractGST(25,2)
GST25B3.df <- subtractGST(25,3)

GST6.25B1.df <- subtractGST(6.25,1)
GST6.25B2.df <- subtractGST(6.25,2)
GST6.25B3.df <- subtractGST(6.25,3)

GST1.56B1.df <- subtractGST(1.56,1)
GST1.56B2.df <- subtractGST(1.56,2)
GST1.56B3.df <- subtractGST(1.56,3)

GST0.39B1.df <- subtractGST(0.39,1)
GST0.39B2.df <- subtractGST(0.39,2)
GST0.39B3.df <- subtractGST(0.39,3)

GST0.10B1.df <- subtractGST(0.10, 1)
GST0.10B2.df <- subtractGST(0.10, 2)
GST0.10B3.df <- subtractGST(0.10, 3)

#combine all data back together again
#filter out the GST tagged antigens
no_tags.df <- filter(target.df, !(Expression_Tag == "GST"))
no_tags.df <- tibble::column_to_rownames(no_tags.df, var="Name")
no_tags.df <- no_tags.df[,(ncol(no_tags.df) - ncol(t(trans.norm.avg))+1):ncol(no_tags.df)]
#then rbind the GST data frames to that one. The order of the targets shouldn't matter anymore
optimization.df <- rbind(no_tags.df, GST0.10B1.df, GST0.10B2.df, GST0.10B3.df,
    GST0.39B1.df, GST0.39B2.df, GST0.39B3.df, GST1.56B1.df, GST1.56B2.df, GST1.56B3.df, 
    GST100B1.df, GST100B2.df, GST100B3.df, GST25B1.df, GST25B2.df, GST25B3.df,
    GST6.25B1.df, GST6.25B2.df, GST6.25B3.df)
#Check that dimensions are the same as before
dim(optimization.df) == dim(t(trans.norm.avg))
View(optimization.df)

# ratio of positive to negative - subtract the negative log2 value for each condition
optimization.df <- t(optimization.df)

Neg <- optimization.df[c(grep("Neg", rownames(optimization.df))),]
CP3 <- optimization.df[c(grep("CP3", rownames(optimization.df))),]
PRISM <- optimization.df[c(grep("PRISM", rownames(optimization.df))),]
Swazi <- optimization.df[c(grep("Swazi", rownames(optimization.df))),]

subtractNeg <- function(x){
for(i in 1:nrow(Neg)){
  for(k in 1:ncol(Neg)){
    if(is.na(Neg[i,k])|is.na(x[i,k])){
      x[i,k] <- NA
    } else {
      x[i,k] <- x[i,k] - Neg[i,k]
    }
  }
}
return(x)
}

#Now a negative value means the positive is less than the negative control
#Do not set these values to 0
CP3neg <- subtractNeg(CP3)
PRISMneg <- subtractNeg(PRISM)
Swazineg <- subtractNeg(Swazi)

Finaldata <- rbind(CP3neg,PRISMneg,Swazineg)

#get the data in the right format for ggplot2 and linear mixed models
#this doesn't include the print buffer as a factor
finaldata.df <- merge(sample_meta_f.df, Finaldata, by.y = "row.names", by.x = "sample_id", sort = TRUE)
#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(finaldata.df))
colnames(finaldata.df) <- newnames

#how to get print buffer as a factor for each antigen
#merge transposed final data with target metadata again
#separate data by print buffer in three data frames
#replace row names with new rownames that are the same for all print buffers, don't include spaces
#bind these back together.

#at some point need to change the sample_id to only CP3, PRISM, and Swazi 
#so that the model knows they are the same sample

#Print buffers: 1 = AJ Buffer C		2 = AJ Glycerol buffer		3 = Nexterion Spot

#once the data is formatted right, do linear mixed model for each antigen
#call the data frame optimization
#fixed effects are: slide_type, print buffer (PB), blocking buffer (BB), blocking buffer dilution (BBD)
#random effects are: sample; using random slope model where print buffer and blocking buffer (and dilution) affect sample variation
colnames(finaldata.df[23])
#AMA1, 100 ug/mL, print buffer 1
fullmodel <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + (1+blocking_buffer+block_dilution|sample_id), REML = FALSE, data = finaldata.df)
