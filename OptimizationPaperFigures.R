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

#################################################
######### Figure 5 ##############################
#################################################

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


