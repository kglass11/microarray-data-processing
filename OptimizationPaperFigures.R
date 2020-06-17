#Optimization figures final for paper 

#Consolidate data input and figure generation from the 3 previous scripts
#Only including plots that will go into the paper

#data input
setwd("/Users/Katie/Desktop/R files from work/OptimizationPaperFiguresFinal")
getwd()

library(broom)
library(lme4)
library(multcomp)
library(emmeans)
library(pbkrtest)
library(phia)
library(lmerTest)
library(multcompView)
library(MuMIn)

require("gtools")

library(Rmisc)
library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)
library(NeatMap)

library(corrgram)
library(corrplot)

#################################################
######### Figure 2 ##############################
#################################################

#Need to get plot script from Tate
#This will get moved to supplemental data?

#################################################
#### Figure 3 - Identify Best Conditions ########
#################################################

#drop Figure 3B
#organize top 8 conditions by slide type in alpha order

load("OptimizationLMMready.RData") #this file is generated 
#from the script "optimizationanalysis.R" prior to running the LMM models
#the script which has the original plots for figure 3 is "OptimizationTop10v2.R"
#which just says it continues on from the optimizationanalysis.R script.

#################################################
#### Figure 4 - Replicate Correlations###########
#################################################
rm(list=ls())

#copied from OptimizationRepCorrelation.R
load("RepsCor.RData")

#Correlation replicate 1 with replicate 2 for all conditions (480 conditions) 

#rep1.matrix and rep2.matrix have the data with negative normalized values set to zero
#trans.norm.rep1 and trans.norm.rep2 have the data without setting negative normalized values to zero
#trans.norm.meta.rep1 and 2 have the data with sample_id as a variable not rownames and have 
#slide_type, blocking_buffer, and block_dilution as final 3 columns
#data is not a ratio of positive to negative, nor has GST been subtracted
#data includes buffer and ref spots as well, but exclude these when separating by print buffer.
#We don't want to count those because they are printed in system buffer only, and don't have a specific print buffer
#Separate data by print buffers for replicate 1 and replicate 2.

#make a target_meta.df with the block 2 antigen names
target_meta2.df <- target_meta.df
target_meta2.df$Name <- colnames(trans.norm.rep2)

#merge transposed final data with target metadata again
target1.df <- merge(target_meta.df, t(trans.norm.rep1), by.y="row.names", by.x = "Name", sort = FALSE)
target2.df <- merge(target_meta2.df, t(trans.norm.rep2), by.y="row.names", by.x = "Name", sort = FALSE)

#separate data by print buffer in three data frames
PB1.1 <- filter(target1.df, Print_Buffer == "1")
PB2.1 <- filter(target1.df, Print_Buffer == "2")
PB3.1 <- filter(target1.df, Print_Buffer == "3")

PB1.2 <- filter(target2.df, Print_Buffer == "1")
PB2.2 <- filter(target2.df, Print_Buffer == "2")
PB3.2 <- filter(target2.df, Print_Buffer == "3")

#replace row names with new rownames that are the same for all print buffers, don't include spaces
rownames(PB1.1) <- PB1.1$Name
rownames(PB2.1) <- PB1.1$Name
rownames(PB3.1) <- PB1.1$Name

rownames(PB1.2) <- PB1.2$Name
rownames(PB2.2) <- PB1.2$Name
rownames(PB3.2) <- PB1.2$Name

#Print buffers: 1 = AJ Buffer C		2 = AJ Glycerol buffer		3 = Nexterion Spot
PB1.1 <- as.data.frame(t(PB1.1[,6:ncol(PB1.1)]))
PB1.1$print_buffer <- "AJ"
PB2.1 <- as.data.frame(t(PB2.1[,6:ncol(PB2.1)]))
PB2.1$print_buffer <- "AJ_Glycerol"
PB3.1 <- as.data.frame(t(PB3.1[,6:ncol(PB3.1)]))
PB3.1$print_buffer <- "Nexterion_Spot"

PB1.2 <- as.data.frame(t(PB1.2[,6:ncol(PB1.2)]))
PB1.2$print_buffer <- "AJ"
PB2.2 <- as.data.frame(t(PB2.2[,6:ncol(PB2.2)]))
PB2.2$print_buffer <- "AJ_Glycerol"
PB3.2 <- as.data.frame(t(PB3.2[,6:ncol(PB3.2)]))
PB3.2$print_buffer <- "Nexterion_Spot"

#merge these with sample_meta_f to get other columns, then rbind all together
PB1.1.meta <- merge(sample_meta_f.df, PB1.1, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB2.1.meta <- merge(sample_meta_f.df, PB2.1, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB3.1.meta <- merge(sample_meta_f.df, PB3.1, by.x = "sample_id", by.y = "row.names",sort=FALSE)

PB1.2.meta <- merge(sample_meta_f.df, PB1.2, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB2.2.meta <- merge(sample_meta_f.df, PB2.2, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB3.2.meta <- merge(sample_meta_f.df, PB3.2, by.x = "sample_id", by.y = "row.names",sort=FALSE)

Rep1.df <- rbind(PB1.1.meta, PB2.1.meta, PB3.1.meta)
Rep2.df <- rbind(PB1.2.meta, PB2.2.meta, PB3.2.meta)

#need to change the sample_id to only CP3, PRISM, and Swazi 
#so that the model knows they are the same sample
Rep1.df$sample <- "CP3"
Rep1.df$sample[c(grep("PRISM", Rep1.df$sample_id))] <- "PRISM"
Rep1.df$sample[c(grep("Swazi", Rep1.df$sample_id))] <- "Swazi"
Rep1.df$sample[c(grep("Neg", Rep1.df$sample_id))] <- "Neg"

Rep2.df$sample <- "CP3"
Rep2.df$sample[c(grep("PRISM", Rep2.df$sample_id))] <- "PRISM"
Rep2.df$sample[c(grep("Swazi", Rep2.df$sample_id))] <- "Swazi"
Rep2.df$sample[c(grep("Neg", Rep2.df$sample_id))] <- "Neg"

#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(Rep1.df))
colnames(Rep1.df) <- newnames

newnames <- make.names(colnames(Rep2.df))
colnames(Rep2.df) <- newnames

#calculate pearson's correlation coefficients for replicate 1 vs replicate 2 row-wise
#there are 30 rows for which all columns have NA, so correlation matrix also has NA for those conditions
repcor <- c(nrow(Rep1.df))

for(i in 1:nrow(Rep1.df)){
  repcor[i] <- cor(x = c(as.matrix(Rep1.df[i,11:46])), y = c(as.matrix(Rep2.df[i,11:46])), use = "everything")
}
remove(i)

#add correlation coefficients to the metadata, but have correlation coeffs for Neg data as well
repcor.df <- as.data.frame(repcor)
repcor.df$sample_id <- Rep1.df$sample_id
repcor.meta.df <- merge(sample_meta_f.df, repcor.df, sort = TRUE)

#add print buffer and sample - *** need to make sure these are the right print buffer and sample!!
sortRep1 <- Rep1.df[order(Rep1.df$sample_id),]
repcor.meta.df$print_buffer <- sortRep1$print_buffer
repcor.meta.df$sample <- sortRep1$sample

#set factor order for sample
repcor.meta.df$sample <- factor(repcor.meta.df$sample, levels = as.character(c("CP3", "PRISM", "Swazi", "Neg")))

#CP3alone4 has the 8 common Top 10% conditions
rep8 <- merge(CP3alone4, repcor.meta.df)
rep8melt <- melt(rep8)
rep8sub <- filter(rep8melt, variable == "repcor")

rep8sub$rowid <- factor(rep8sub$rowid, levels = rev(unique(Top8rowidorder)))

#Updated Figure 4A. Plot correlation coefficients vs. slide type, CP3 ONLY + points for Top 8 Conditions
repcor.CP3 <- filter(repcor.meta.df, sample == "CP3")
rep8CP3 <- filter(rep8sub, sample == "CP3")

#set the seed so that the black dots will be in the same spot every time
#seed 9 is ok but can look for a better one later if we want to
set.seed(9)

png(filename = paste0("Fig4A.CP3repcor.tif"), width = 7.5, height = 3.5, units = "in", res = 1200)

ggplot(data = repcor.CP3, aes(x = slide_type, y=repcor)) + geom_jitter(size = 0.6, height = 0, width = 0.33, aes(color = slide_type)) + 
  geom_jitter(data = rep8CP3, aes(x=slide_type, y = value), size = 0.6, width = 0.45, height = 0, color = "black") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) + 
  theme(legend.position = "none") +
  labs(x = "Slide Type", y = "Pearson's Correlation Coefficient")

graphics.off() 

#Figure 4B
#CP3 "Best" conditions for 4 antigens only - correlogram
#data frame hello has all the right data, but may need to transpose and subset

hellocor <- as.data.frame(t(hello[,13:48]))
colnames(hellocor) <- hello$rowid

#reorder the columns in hellocor to match Top8rowidorder
corOrd <- hellocor[,Top8rowidorder]

#plot with corrplot package
png(filename = paste0("Fig4B.Common8.CorrelogramV2.tif"), width = 5, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

corrplot.mixed(cor(corOrd, use = "complete.obs"), tl.col="black")

graphics.off()

#calculate how many arrays have correlation coefficients above certain thresholds
#function to calculate the percent of correlation coefficients above a certain threshold
repcor.per <- function(threshold, vector){
  pepper <- (vector > threshold) + 0
  x <- sum(pepper, na.rm=TRUE)/length(vector)*100
  return(paste(round(x,1), "% of correlations coefficients are above", threshold))
}

repcor.per(0.95, repcor)
repcor.per(0.9, repcor)
repcor.per(0.8, repcor)
repcor.per(0.7, repcor)
repcor.per(0.5, repcor)

#repeat calculations but for CP3 only
CP3repcor.df <- filter(repcor.meta.df, sample == "CP3")
CP3repcor <- c(CP3repcor.df$repcor)

repcor.per(0.95, CP3repcor)
repcor.per(0.9, CP3repcor)
repcor.per(0.8, CP3repcor)
repcor.per(0.7, CP3repcor)
repcor.per(0.5, CP3repcor)



#################################################
#### Figure 5 - Repeats Dilution Curves #########
#################################################
rm(list=ls())

#copied from optimization analysis repeats script (optimizationanalysisREPEATS.R)
load("OptRepeats_AfterProcessing.RData")

#Data processing for figure 5

#ratio of positive to negative - subtract the negative log2 value for each condition
optimization.df <- t(optimization.df)

#make individual data frames for each sample
Neg <- optimization.df[c(grep("Neg", rownames(optimization.df))),]
CP3 <- optimization.df[c(grep("CP3", rownames(optimization.df))),]
NIBSC <- optimization.df[c(grep("10/198", rownames(optimization.df))),]

#Function for subtracting negatives to get ratio 
#Note: there is one sample missing from NIBSC - but it's the last sample (#24), so the function should still work correctly 
subtractNeg <- function(x){
  for(i in 1:nrow(x)){
    for(k in 1:ncol(x)){
      if(is.na(Neg[i,k])|is.na(x[i,k])){
        x[i,k] <- NA
      } else {
        x[i,k] <- x[i,k] - Neg[i,k]
      }
    }
  }
  return(x)
}


#GST values for each sample
CP3.GST <- CP3[,c(grep("GST", colnames(CP3)))]
NIBSC.GST <- NIBSC[,c(grep("GST", colnames(NIBSC)))]
Neg.GST <- Neg[,c(grep("GST", colnames(Neg)))]

#Now a negative value means the positive is less than the negative control
#Do not set these values to 0
CP3neg <- subtractNeg(CP3)
NIBSCneg <- subtractNeg(NIBSC)

Finaldata <- rbind(CP3neg,NIBSCneg)

#get the data in the right format for ggplot2 (also for linear mixed models, but we aren't doing this right now)

#how to get print buffer as a factor for each antigen:
#merge transposed final data with target metadata again
#separate data by print buffer in three data frames
#replace row names with new rownames that are the same for all print buffers, don't include spaces
#bind these back together.
finaltarget.df <- merge(target_meta.df, t(Finaldata), by.y="row.names", by.x = "Name", sort = FALSE)
PB1 <- filter(finaltarget.df, Print_Buffer == "1")
PB2 <- filter(finaltarget.df, Print_Buffer == "2")
PB3 <- filter(finaltarget.df, Print_Buffer == "3")

rownames(PB1) <- PB1$Name
rownames(PB2) <- PB1$Name
rownames(PB3) <- PB1$Name

#Print buffers: 1 = AJ Buffer C		2 = AJ Glycerol buffer		3 = Nexterion Spot
PB1 <- as.data.frame(t(PB1[,6:ncol(PB1)]))
PB1$print_buffer <- "AJ"
PB2 <- as.data.frame(t(PB2[,6:ncol(PB2)]))
PB2$print_buffer <- "AJ_Glycerol"
PB3 <- as.data.frame(t(PB3[,6:ncol(PB3)]))
PB3$print_buffer <- "Nexterion_Spot"

#merge these with sample_meta_f to get other columns, then rbind all together
PB1_meta <- merge(sample_meta_f.df, PB1, by.x = "sample_id_unique", by.y = "row.names",sort=FALSE)
PB2_meta <- merge(sample_meta_f.df, PB2, by.x = "sample_id_unique", by.y = "row.names",sort=FALSE)
PB3_meta <- merge(sample_meta_f.df, PB3, by.x = "sample_id_unique", by.y = "row.names",sort=FALSE)

AHHHH.df <- rbind(PB1_meta, PB2_meta, PB3_meta)

#need to add a sample identifier of only CP3 and NIBSC (10/198)
#so that the model knows they are the same sample
AHHHH.df$sample <- "CP3"
AHHHH.df$sample[c(grep("10/198", AHHHH.df$sample_id))] <- "NIBSC"

#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(AHHHH.df))
colnames(AHHHH.df) <- newnames

#Make subsets which are based on CP3 or NIBSC alone
AHHsub <- filter(AHHHH.df, sample == "CP3" | sample == "NIBSC")
CP3all <- filter(AHHHH.df, sample == "CP3")
CP3sub <- CP3all[12:47]
NIBSCall <- filter(AHHHH.df, sample == "NIBSC")
NIBSCsub <- NIBSCall[12:47]

AHHsub2 <- rbind(NIBSCall, CP3all)

conditions <- CP3all[,c(9,10,11,48)]
#isolate the 18 unique conditions
cond.unique <- unique(conditions)
#get an identifying number for each condition
cond.unique <- tibble::rowid_to_column(cond.unique)

#export conditions and identifying numbers
write.csv(cond.unique,file = "RowIDKey.csv")

#get a data frame with all of the data and the cond.unique info + the sample column
AHHsub3 <- merge(cond.unique, AHHsub2[,c(2,9:49)], sort = FALSE)
AHHsub3$rowid <- as.character(AHHsub3$rowid)

#melt and subset data for the plots
melt1 <- melt(AHHsub3, variable.name = "Antigen")

#Dilution curve plots with error bars! 
  
#Prepare data frame with concentration information
meltcon <- melt1
meltcon$concentration <- 100

meltcon$concentration[c(grep("25", meltcon$Antigen))] <- 25
meltcon$concentration[c(grep("6.25", meltcon$Antigen))] <- 6.25
meltcon$concentration[c(grep("1.56", meltcon$Antigen))] <- 1.56
meltcon$concentration[c(grep("0.39", meltcon$Antigen))] <- 0.39
meltcon$concentration[c(grep("0.1", meltcon$Antigen))] <- 0.1

#Prepare summary data frame for error bars
meltconsum <- summarySE(meltcon, na.rm=TRUE, measurevar = "value", groupvars = c("Antigen", "sample", "rowid", "concentration"))

#MSP1-19
meltMSP <- meltconsum[c(grep("MSP1", meltconsum$Antigen)),]
meltMSP$rowid <- factor(meltMSP$rowid, levels = 1:18)

#AMA1
meltAMA1 <- meltconsum[c(grep("AMA1", meltconsum$Antigen)),]
meltAMA1$rowid <- factor(meltAMA1$rowid, levels = 1:18)

#paper figure 5 - AMA1 and MSP1 dilutions for the previously selected 8 conditions
#condition numbers from this study: 2,3,4,5,6,7,13,15
#Later, am going to have rename the conditions for all the figures to 1-8, once we decide what order we want them in

meltMSPsub2 <- filter(meltMSP, rowid %in% c(2,3,4,5,6,7,13,15))

png(filename = "NIBSC.MSP.Dil.Fig5.tif", width = 4, height = 3.3, units = "in", res = 1200)

print(ggplot(filter(meltMSPsub2, sample == "NIBSC"), aes(x = as.factor(concentration), y = as.numeric(value), color = rowid)) + 
        geom_errorbar(aes(ymin=value-se, ymax=value+se), color="black", width=.1) +
        geom_point(shape=18, size = 2) +
        geom_line(aes(group = rowid)) + 
        theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Concentration (µg/mL)") + 
        theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        scale_color_hue(name = "Condition"))

graphics.off()

png(filename = "CP3.MSP.Dil.Fig5.tif", width = 4, height = 3.3, units = "in", res = 1200)

print(ggplot(filter(meltMSPsub2, sample == "CP3"), aes(x = as.factor(concentration), y = as.numeric(value), color = rowid)) + 
        geom_errorbar(aes(ymin=value-se, ymax=value+se), color="black", width=.1) +
        geom_point(shape=18, size = 2) +
        geom_line(aes(group = rowid)) + 
        theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Concentration (µg/mL)") + 
        theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        scale_color_hue(name = "Condition"))

graphics.off()

meltAMA1sub2 <- filter(meltAMA1, rowid %in% c(2,3,4,5,6,7,13,15))

png(filename = "NIBSC.AMA1.DilFig5.tif", width = 4, height = 3.3, units = "in", res = 1200)

print(ggplot(filter(meltAMA1sub2, sample == "NIBSC"), aes(x = as.factor(concentration), y = as.numeric(value), color = rowid)) + 
        geom_errorbar(aes(ymin=value-se, ymax=value+se), color="black", width=.1) +
        geom_point(shape=18, size = 2) +
        geom_line(aes(group = rowid)) + 
        theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Concentration (µg/mL)") + 
        theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        scale_color_hue(name = "Condition"))

graphics.off()

png(filename = "CP3.AMA1.DilFig5.tif", width = 4, height = 3.3, units = "in", res = 1200)

print(ggplot(filter(meltAMA1sub2, sample == "CP3"), aes(x = as.factor(concentration), y = as.numeric(value), color = rowid)) + 
        geom_errorbar(aes(ymin=value-se, ymax=value+se), color="black", width=.1) +
        geom_point(shape=18, size = 2) +
        geom_line(aes(group = rowid)) + 
        theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Concentration (µg/mL)") + 
        theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        scale_color_hue(name = "Condition"))

graphics.off()


