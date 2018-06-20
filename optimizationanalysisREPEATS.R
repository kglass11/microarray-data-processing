#optimization data analysis - REPEATS!!
#KG
#June 19, 2018

#"I:/Drakeley Group/Protein microarrays/Experiments/270717 Optimisation"
#/Users/Katie/Desktop/R files from work/270717 Optimisation
setwd("/Users/Katie/Desktop/R files from work/170518 Optimisation CP3 repeats")
getwd()

# install.packages("lme4")
# install.packages("broom")
# install.packages("multcomp")
# install.packages("lsmeans")
# install.packages("emmeans")
# install.packages("pbkrtest")
# install.packages("phia")
# install.packages("lmerTest")
# install.packages("multcompView")
# install.packages("MuMIn")
# install.packages("NeatMap")

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

load("OptRepeats_AfterProcessing.RData")

# ratio of positive to negative - subtract the negative log2 value for each condition
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


#Plot GST values for each sample
CP3.GST <- CP3[,c(grep("GST", colnames(CP3)))]
NIBSC.GST <- NIBSC[,c(grep("GST", colnames(NIBSC)))]
Neg.GST <- Neg[,c(grep("GST", colnames(Neg)))]

png(filename = paste0(study, "_GST_Sample.tif"), width = 7.5, height = 9, units = "in", res = 600)
par(mfrow=c(2,2), oma=c(4,1,1,1),mar=c(6.1,4.1,1.1,1.1))
boxplot(CP3.GST, pch='*', col = "light blue", ylim=c(0,1.5), main = "CP3 GST",
     ylab="Normalized log2(MFI)", las=2, cex.axis = 0.5)

boxplot(NIBSC.GST, pch='*', col = "light blue", ylim=c(0,1.5), main = "NIBSC GST",
     ylab="Normalized log2(MFI)", las=2, cex.axis = 0.5)

boxplot(Neg.GST, pch='*', col = "light blue", ylim=c(0,1.5), main = "Neg GST",
     ylab="Normalized log2(MFI)", las=2, cex.axis = 0.5)

#print text on the plots for number of samples where GST is above buffer (0)
mtext(paste("GST > 0: CP3 =", round(sum(CP3.GST > 0, na.rm = TRUE), digits=2), 
    "(", round(sum(CP3.GST > 0, na.rm = TRUE)/length(CP3.GST)*100, digits=2), "%), NIBSC = ",
    round(sum(NIBSC.GST > 0, na.rm = TRUE), digits=2), 
    "(", round(sum(NIBSC.GST > 0, na.rm = TRUE)/length(NIBSC.GST)*100, digits=2), "%), Neg = ",
    round(sum(Neg.GST > 0, na.rm = TRUE), digits=2),
    "(", round(sum(Neg.GST > 0, na.rm = TRUE)/length(Neg.GST)*100, digits=2), "%)"), 
    side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Now a negative value means the positive is less than the negative control
#Do not set these values to 0
CP3neg <- subtractNeg(CP3)
NIBSCneg <- subtractNeg(NIBSC)

Finaldata <- rbind(CP3neg,NIBSCneg)

#get the data in the right format for ggplot2 and linear mixed models

#how to get print buffer as a factor for each antigen
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
PB1_meta <- merge(sample_meta_f.df, PB1, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB2_meta <- merge(sample_meta_f.df, PB2, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB3_meta <- merge(sample_meta_f.df, PB3, by.x = "sample_id", by.y = "row.names",sort=FALSE)

AHHHH.df <- rbind(PB1_meta, PB2_meta, PB3_meta)

#need to add a sample identifier of only CP3, PRISM, and Swazi 
#so that the model knows they are the same sample
AHHHH.df$sample <- "CP3"
AHHHH.df$sample[c(grep("PRISM", AHHHH.df$sample_id))] <- "PRISM"
AHHHH.df$sample[c(grep("Swazi", AHHHH.df$sample_id))] <- "Swazi"

#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(AHHHH.df))
colnames(AHHHH.df) <- newnames

#####Plotting#######

#Plot the data for the PRISM and CP3 samples only, organized by geometric mean

AHHsub <- filter(AHHHH.df, sample == "CP3" | sample == "PRISM")
CP3all <- filter(AHHHH.df, sample == "CP3")
CP3sub <- CP3all[11:46]
PRISMall <- filter(AHHHH.df, sample == "PRISM")
PRISMsub <- PRISMall[11:46]

geomean <- (CP3sub + PRISMsub)/2
geomean <- tibble:: rowid_to_column(geomean)
PRISMall <- tibble::rowid_to_column(PRISMall)
CP3all <- tibble::rowid_to_column(CP3all)
AHHsub2 <- rbind(PRISMall, CP3all)

#plot data from top 48 geomeans (10%) of conditions; AMA1, 100
geomeanAMA1.100 <- geomean[order(geomean$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
rowidorder <- c(as.character(geomeanAMA1.100$rowid))
geomeanAMA1.100$rowid <- factor(geomeanAMA1.100$rowid, levels = unique(rowidorder))

AHHsub2$rowid <- factor(AHHsub2$rowid, levels = unique(rowidorder))

conditions <- PRISMall[,c(1,9,10,11,48)]

#export table of conditions matching top 10% rowid. 
conditionsAMA1.100 <- merge(geomeanAMA1.100, conditions, sort = FALSE)
write.csv(conditionsAMA1.100[1:48, c(1,38:41)], file = "top10AMA1.100.csv")

png(filename = paste0("AMA1.100.GM.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanAMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

geomeanAMA1.100$rowid <- factor(geomeanAMA1.100$rowid, levels = rev(unique(rowidorder)))

png(filename = paste0("V.AMA1.100.GM.tif"), width = 3.5, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanAMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 7)) +
  coord_flip() + labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "AMA1 100 ?g/mL") +
  ylim(0,8)

graphics.off()

  #plot the original CP3 and PRISM values instead of the geomean
  #keep the rowid order from the geomean

  GMcondAMA1 <- conditionsAMA1.100[1:48, c(1,38:41)]
  DataGM.AMA1 <- merge(GMcondAMA1, AHHsub2, sort=FALSE)
  
  png(filename = paste0("AMA1.samples.GM.tif"), width = 8, height = 3.5, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(DataGM.AMA1, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample))  + 
    geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "AMA1 100 µg/mL")
  
  graphics.off()  


#MSP1-19, 100 ug/mL; plot data from top 48 geomeans (10%) of conditions;
colnames(geomean[8])
geomeanMSP1.19.100 <- geomean[order(geomean$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
rowidorder <- c(as.character(geomeanMSP1.19.100$rowid))
geomeanMSP1.19.100$rowid <- factor(geomeanMSP1.19.100$rowid, levels = unique(rowidorder))

AHHsub2$rowid <- factor(AHHsub2$rowid, levels = unique(rowidorder))

conditions <- PRISMall[,c(1,9,10,11,48)]

#export table of conditions matching top 10% rowid. 
conditionsMSP1.19.100 <- merge(geomeanMSP1.19.100, conditions, sort = FALSE)
write.csv(conditionsMSP1.19.100[1:48, c(1,38:41)], file = "top10MSP1.19.100.csv")

png(filename = paste0("MSP1.19.100.GM.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanMSP1.19.100[1:48,], aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

geomeanMSP1.19.100$rowid <- factor(geomeanMSP1.19.100$rowid, levels = rev(unique(rowidorder)))

png(filename = paste0("V.MSP1.19.100.GM.tif"), width = 3.5, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanMSP1.19.100[1:48,], aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 7)) +
  coord_flip() + labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "MSP1-19 100 ?g/mL") +
  ylim(0,8)

graphics.off()

#plot the original CP3 and PRISM values instead of the geomean
#keep the rowid order from the geomean

GMcondMSP1 <- conditionsMSP1.19.100[1:48, c(1,38:41)]
DataGM.MSP1 <- merge(GMcondMSP1, AHHsub2, sort=FALSE)

png(filename = paste0("MSP1.19.samples.GM.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(DataGM.MSP1, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample))  + 
  geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "MSP1-19 100 µg/mL")

graphics.off()  

#Hyp2, 100 ug/mL; plot data from top 48 geomeans (10%) of conditions;
colnames(geomean[14])
geomeanHyp2.100 <- geomean[order(geomean$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
rowidorder <- c(as.character(geomeanHyp2.100$rowid))
geomeanHyp2.100$rowid <- factor(geomeanHyp2.100$rowid, levels = unique(rowidorder))

AHHsub2$rowid <- factor(AHHsub2$rowid, levels = unique(rowidorder))

conditions <- PRISMall[,c(1,9,10,11,48)]

#export table of conditions matching top 10% rowid. 
conditionsHyp2.100 <- merge(geomeanHyp2.100, conditions, sort = FALSE)
write.csv(conditionsHyp2.100[1:48, c(1,38:41)], file = "top10Hyp2.100.csv")

png(filename = paste0("Hyp2.100.GM.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanHyp2.100[1:48,], aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

geomeanHyp2.100$rowid <- factor(geomeanHyp2.100$rowid, levels = rev(unique(rowidorder)))

png(filename = paste0("V.Hyp2.100.GM.tif"), width = 3.5, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanHyp2.100[1:48,], aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 7)) +
  coord_flip() + labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Hyp2 100 ?g/mL") +
  ylim(0,8)

graphics.off()

#plot the original CP3 and PRISM values instead of the geomean
#keep the rowid order from the geomean

GMcondHyp2 <- conditionsHyp2.100[1:48, c(1,38:41)]
DataGM.Hyp2 <- merge(GMcondHyp2, AHHsub2, sort=FALSE)

png(filename = paste0("Hyp2.samples.GM.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(DataGM.Hyp2, aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1, color = sample))  + 
  geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Hyp2 100 µg/mL")

graphics.off()  

#EPF1v2, 100 ug/mL; plot data from top 48 geomeans (10%) of conditions;
colnames(geomean[26])
geomeanEPF1v2.100 <- geomean[order(geomean$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
rowidorder <- c(as.character(geomeanEPF1v2.100$rowid))
geomeanEPF1v2.100$rowid <- factor(geomeanEPF1v2.100$rowid, levels = unique(rowidorder))

AHHsub2$rowid <- factor(AHHsub2$rowid, levels = unique(rowidorder))

conditions <- PRISMall[,c(1,9,10,11,48)]

#export table of conditions matching top 10% rowid. 
conditionsEPF1v2.100 <- merge(geomeanEPF1v2.100, conditions, sort = FALSE)
write.csv(conditionsEPF1v2.100[1:48, c(1,38:41)], file = "top10EPF1v2.100.csv")

png(filename = paste0("EPF1v2.100.GM.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanEPF1v2.100[1:48,], aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

geomeanEPF1v2.100$rowid <- factor(geomeanEPF1v2.100$rowid, levels = rev(unique(rowidorder)))

png(filename = paste0("V.EPF1v2.100.GM.tif"), width = 3.5, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanEPF1v2.100[1:48,], aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 7)) +
  coord_flip() + labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "EPF1v2 100 ?g/mL") +
  ylim(0,8)

graphics.off()

GMcondEPF1v2 <- conditionsEPF1v2.100[1:48, c(1,38:41)]
DataGM.EPF1v2 <- merge(GMcondEPF1v2, AHHsub2, sort=FALSE)

png(filename = paste0("EPF1v2.samples.GM.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(DataGM.EPF1v2, aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1, color = sample))  + 
  geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "EPF1v2 100 µg/mL")

graphics.off()  

#determine which conditions are in the top 10% for AMA1, MSP1-19, Hyp2, and EPF1v2

conditionAMA1sub <- conditionsAMA1.100[1:48, c(1,38:41)]
conditionMSP1.19sub <- conditionsMSP1.19.100[1:48, c(1,38:41)]
conditionHyp2sub <- conditionsHyp2.100[1:48, c(1,38:41)]
conditionEPF1v2sub <- conditionsEPF1v2.100[1:48, c(1,38:41)]

two <- merge(conditionAMA1sub, conditionMSP1.19sub, sort = FALSE)
three <- merge(two, conditionHyp2sub, sort = FALSE)
four <- merge(three, conditionEPF1v2sub, sort = FALSE)

write.csv(four, file = "commontop10.csv")

#plot all geometric means for all 4 antigens for common 11 conditions
#currently row ID are sorted based on the AMA1 geomeans

fourGM <- merge(four, geomean, sort = FALSE)
fourSamples <- merge(four, AHHsub2, sort = FALSE)

fourmelt <- melt(fourGM)
foursamplesmelt <- melt(fourSamples)

fourmeltsub <- filter(fourmelt, variable == "X37_1.EPF1v2.100ug.ml_1" | variable == "X25_1.Hyp2.100ug.ml_1"
                      | variable == "X19_1.PfMSP1.19.100ug.ml_1" | variable == "X13_1.PfAMA1.100ug.ml_1" )

foursamplessub <- filter(foursamplesmelt, variable == "X37_1.EPF1v2.100ug.ml_1" | variable == "X25_1.Hyp2.100ug.ml_1"
                         | variable == "X19_1.PfMSP1.19.100ug.ml_1" | variable == "X13_1.PfAMA1.100ug.ml_1" )

png(filename = paste0("GMfourAgstop10.tif"), width = 5, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(fourmeltsub, aes(x = rowid, y = value, color = variable)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 10)) +
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Geometric mean of common combinations") +
  ylim(0,4) + scale_color_hue(labels = c("AMA1", "MSP1-19", "Hyp2", "EPF1v2"))

graphics.off()

#vertical plot of the same thing
rowids <- c(as.character(fourGM$rowid))
fourmeltsub$rowid <- factor(fourmeltsub$rowid, levels = rev(unique(rowids)))
foursamplessub$rowid <- factor(foursamplessub$rowid, levels = rev(unique(rowids)))

png(filename = paste0("V.GMfourAgstop10.tif"), width = 5, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(fourmeltsub, aes(x = rowid, y = value, color = variable)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 10)) +
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Geometric mean of common combinations") +
  ylim(0,4) + scale_color_hue(labels = c("AMA1", "MSP1-19", "Hyp2", "EPF1v2")) + coord_flip()

graphics.off()

#Plot the values for CP3 and PRISM for the top combinations rather than geomean
png(filename = paste0("V.samples.top10.tif"), width = 5, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(foursamplessub, aes(x = rowid, y = value, color = variable, shape = sample)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 10)) +
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 and PRISM values for common combinations") +
  scale_color_hue(labels = c("AMA1", "MSP1-19", "Hyp2", "EPF1v2")) + coord_flip()

graphics.off()



#plotting samples instead of geomean

png(filename = paste0("AMA1.100.2samples.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHsub2, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()




#violin, scatter, or boxplots of the data

#all samples and all variables - AMA1 100 - Scatter plot
png(filename = paste0(study, "AMA1.100_facet_grid.tif"), width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHHH.df, aes(x=slide_type, y=X13_1.PfAMA1.100ug.ml_1, color = sample, shape = blocking_buffer)) + 
  geom_point(size = 1)  + 
  theme_bw() + 
  theme(text = element_text(size=9), axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.1)) + 
  labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = "AMA1 (100 µg/mL)") +
  facet_grid(print_buffer ~ block_dilution )

graphics.off()

#violin plot of slide type and print buffer, everything else combined 
png(filename = paste0(study, "AMA1.100_ST_PB.tif"), width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHHH.df, aes(x=reorder(slide_type, X13_1.PfAMA1.100ug.ml_1, FUN=max), y=X13_1.PfAMA1.100ug.ml_1, color = print_buffer, fill = print_buffer)) + 
  geom_violin(trim = FALSE) + 
  theme_bw() +
  guides(color=FALSE) +
  theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, hjust = 1), legend.position="top") +
  labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = "AMA1 (100 µg/mL)", fill = "Print Buffer")

graphics.off()

#Make a function of the plots to use with other columns

#column is the column name from AHHHH.df for each antigen and dilution
#antigen and dilution are character vectors 
# antigenplots <- function(column, antigen, dilution){
#   
# #all samples and all variables - Scatter plot faceted
# png(filename = paste0(antigen,"_", dilution, "_facet_grid.tif"), width = 8, height = 5, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# p1 <- ggplot(AHHHH.df, aes(x=slide_type, y=column, color = sample, shape = blocking_buffer)) + 
#   geom_point(size = 1)  + 
#   theme_bw() + 
#   theme(text = element_text(size=9), axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.1)) + 
#   labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = paste(antigen, dilution, "µg/mL")) +
#   facet_grid(print_buffer ~ block_dilution )
# 
# print(p1)
# graphics.off()
# 
# #violin plot of slide type and print buffer, everything else combined 
# png(filename = paste0(antigen,"_", dilution, "_ST_PB.tif"), width = 8, height = 5, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# p2 <- ggplot(AHHHH.df, aes(x=slide_type, y=column, color = print_buffer, fill = print_buffer)) + 
#   geom_violin(trim = FALSE) + 
#   theme_bw() +
#   guides(color=FALSE) +
#   theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, hjust = 1), legend.position="top") +
#   labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = paste(antigen, dilution, "µg/mL"), fill = "Print Buffer")
# 
# print(p2)
# graphics.off()
# 
# }
# 
# antigenplots(X13_1.PfAMA1.100ug.ml_1, "AMA1", "100")

#heatmap of all the data 
#want a heatmap with antigen name and dilution on the top (x)
#and the combo on the side. Sorted by highest total intensity (sum of the rows)

#prepare relevant data as a matrix
heatdata <- as.matrix(AHHHH.df[,11:(ncol(AHHHH.df)-2)])

#prepare vector for order of rows and columns
#AHHHH.df is currently ordered by sample > print buffer > slide type > blocking buffer > dilution
#and the columns are ordered by antigen > dilution

#prepare character vector labels? or make labels later?
png(filename = paste0(study, "heatmap2.All.tif"), width = 4.5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

heatmap1(heatdata, row.order = NULL)

graphics.off()

#heatmap of the data separated by print buffer 
AJ <- filter(AHHHH.df, print_buffer == "AJ")
AJ.mat <- as.matrix(AJ[,11:46])

AJGly <- filter(AHHHH.df, print_buffer == "AJ_Glycerol")
AJGly.mat <- as.matrix(AJGly[,11:46])

Next <- filter(AHHHH.df, print_buffer == "Nexterion_Spot")
Next.mat <- as.matrix(Next[,11:46])

png(filename = paste0(study, "AJ_heatmap.All.tif"), width = 4.5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
heatmap1(AJ.mat, row.order = NULL)
graphics.off()



  
# 
# heatmap(as.matrix(data))
# 
# #heatmap with ggplot2 - not done yet
# ggplot(AHHHH.df, aes(slide_type, blocking_buffer, block_dilution)) +
#   geom_tile(aes(fill = AHHHH.df[11]), color = "white") +
#   scale_fill_gradient(low = "red", high = "steelblue") +
#   ylab("") +
#   xlab("") +
#   theme(legend.title = element_text(size = 9),
#         legend.text = element_text(size = 9),
#         plot.title = element_text(size=9),
#         axis.title=element_text(size=9,face="bold"),
#         axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(fill = "Normalized Log2(MFI)")
# 
# #hierarchical clustering - this isn't working
# clusters <- hclust(dist(AHHHH.df[11]), na.rm = TRUE)
# plot(clusters)

######### Linear Mixed Effects Models ########

#once the data is formatted right, do linear mixed model for each antigen
#fixed effects are: slide_type, print_buffer, blocking_buffer, block_dilution
#random effects are: sample; using random slope model where print buffer and blocking buffer (and dilution) affect sample variation

#save workspace image to start from this point
save.image(file = "OptimizationLMMready.RData")

#AMA1, 100 ug/mL, - change REML to FALSE for anova
colnames(AHHHH.df[11])
fullmodel <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer 
    + block_dilution + print_buffer + (1|sample), 
    REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#look at the impact of each individual factor by removing 1 factor at a time
#print buffer
model_PB <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer 
               + block_dilution + (1|sample), 
               REML = FALSE, data = AHHHH.df)
  summary(model_PB)
  r.squaredGLMM(model_PB)
  #compare these with likelihood ratio test
  anova(model_PB, fullmodel)
  
#slide type 
  model_ST <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ print_buffer + blocking_buffer 
                   + block_dilution + (1|sample), 
                   REML = FALSE, data = AHHHH.df)
  summary(model_ST)
  r.squaredGLMM(model_ST)
  #compare these with likelihood ratio test
  anova(model_ST, fullmodel)  
  
#blocking buffer
  model_BB <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ print_buffer + slide_type 
                   + block_dilution + (1|sample), 
                   REML = FALSE, data = AHHHH.df)
  summary(model_BB)
  r.squaredGLMM(model_BB)
  #compare these with likelihood ratio test
  anova(model_BB, fullmodel)
  
#block dilution
  model_BD <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ print_buffer + slide_type 
                  + blocking_buffer + (1|sample), 
                   REML = FALSE, data = AHHHH.df)
  summary(model_BD)
  r.squaredGLMM(model_BD)
  #compare these with likelihood ratio test
  anova(model_BD, fullmodel)
  
#look at the 2 way interactions between factors 
  int2model <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                      (slide_type + blocking_buffer + block_dilution + print_buffer)^2 +
                      (1|sample), REML = FALSE, data = AHHHH.df)
  summary(int2model)
  r.squaredGLMM(int2model)
  
  anova(int2model)
  
  #compare these with likelihood ratio test
  anova(fullmodel, int2model)
  
  #Plot residuals - int2model
  png(filename = "AMA1.100.Int2.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int2model),residuals(int2model),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int2model), xlab = "Residuals")
  
  graphics.off()
  
#look at 3 way interactions between factors
  int3model <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                      (slide_type + blocking_buffer + block_dilution + print_buffer)^3 +
                      (1|sample), REML = FALSE, data = AHHHH.df)
  summary(int3model)
  r.squaredGLMM(int3model)
  
  anova(int3model)
  
  #compare these with likelihood ratio test
  anova(fullmodel, int3model)
  anova(int2model, int3model)
  
  #Plot residuals - int3model
  png(filename = "AMA1.100.Int3.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int3model),residuals(int3model),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int3model), xlab = "Residuals")
  
  graphics.off()
  
#look at 4 way interactions between factors
  int4model <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                      (slide_type + blocking_buffer + block_dilution + print_buffer)^4 +
                      (1|sample), REML = FALSE, data = AHHHH.df)
  summary(int4model)
  r.squaredGLMM(int4model)
  
  anova(int4model)
  
  #compare these with likelihood ratio test
  anova(fullmodel, int4model)
  anova(int2model, int4model)
  anova(int3model, int4model)
  
  #Plot residuals - int4model
  png(filename = "AMA1.100.Int4.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int4model),residuals(int4model),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int4model), xlab = "Residuals")
  
  graphics.off()
  
#modified 3 way interaction model - significant 3 way interactions taken from int3model
  int3.2 <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                   slide_type*block_dilution*print_buffer + slide_type*blocking_buffer*block_dilution  + 
                   (1|sample), REML = FALSE, data = AHHHH.df)

  summary(int3.2)
  r.squaredGLMM(int3.2)
  anova(int3.2)
  
  #compare these with likelihood ratio test
  anova(int2model, int3.2)
  anova(int3model, int3.2)
  
  #Plot residuals - int3.2
  png(filename = "AMA1.100.Int3.2.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int3.2),residuals(int3.2),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int3.2), xlab = "Residuals")
  
  graphics.off()
  
  
#prepare a better organized summary of the model and export to a file.
#summary can only be written to a text file, and doesn't keep columns organized
tidy(fullmodel)
augment(fullmodel)

write.csv(tidy(fullmodel), file = "AMA1.100.LMERTidy.csv")
write.csv(augment(fullmodel), file = "AMA1.100.LMERAug.csv")

#test for main effects and interactions with Likelihood Ratio Test?

intmodel <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type * blocking_buffer 
                    * block_dilution * print_buffer + (1|sample), 
                    REML = FALSE, data = AHHHH.df)
summary(intmodel)

#Plot residuals
png(filename = "AMA1.100.Residuals.tif", width = 8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

plot(fitted(intmodel),residuals(intmodel),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
abline(a=0, b=0)
hist(residuals(intmodel), xlab = "Residuals")

graphics.off()


r.squaredGLMM(intmodel)

write.csv(tidy(intmodel), file = "AMA1.100.LMERTidyINT.csv")
write.csv(augment(intmodel), file = "AMA1.100.LMERAugINT.csv")

anova(fullmodel, intmodel)

anova(intmodel)

#get significance and do pairwise comparisons using emmeans
#default adjustment for multiple comparisons is Tukey
#main effects first, however we are not interested in main effects :(
emmeans(fullmodel, pairwise ~ print_buffer)
emmeans(fullmodel, pairwise ~ slide_type)
emmeans(fullmodel, pairwise ~ blocking_buffer)
emmeans(fullmodel, pairwise ~ block_dilution)

#calculate emms
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#repeat for int3.2
emm3.2 <- emmeans(int3.2, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#plot the data
png(filename = "AMA1.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type) + element_text(size=8)

graphics.off()

#plot emmeans for int3.2 model
png(filename = "AMA1.100.emmip3.2.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(emm3.2, print_buffer ~ block_dilution | blocking_buffer | slide_type)

graphics.off()

#extract the emmeans and sort highest to lowest
emm3.2.df <- as.data.frame(emm3.2)
emm3.2.df <- emm3.2.df[order(emm3.2.df$emmean, decreasing = TRUE),]

#make rowid a column to match with AHHHH.df
emm3.2.df <- tibble::rowid_to_column(emm3.2.df)

#merge with AHHHH.df
emmAMA1.100.df <- merge(emm3.2.df, AHHHH.df, all.x = TRUE, sort = FALSE)

#keep same order that we have in the excel file (sorted highest to lowest EMM)
#maybe try to use rownumber as factor, sorted by EMM? and then put in better labels later?
rowidorder <- c(as.character(emmAMA1.100.df$rowid))
emmAMA1.100.df$rowid <- factor(emmAMA1.100.df$rowid, levels = rev(unique(rowidorder)))

#plot original data from groups with higest emmeans (top 36)
emmAMA1.100.df$rowid <- factor(emmAMA1.100.df$rowid, levels = unique(rowidorder))

emmAMA1.100sub <- emmAMA1.100.df[1:108,]

png(filename = paste0("H.AMA1.100.3.2.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(emmAMA1.100sub, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()




#get letters for all pairwise comparisons, interaction model

#this took several hours (>4) to run, and eventually finished. there are too many 
#combinations (480 pairwise), but it worked!!! 
#the contrast function used by cld function automatically uses Tukey adjustment
#for multiple comparisons
pairwiseletters <- cld(emm.intmodel, by = NULL, Letters = LETTERS, sort = TRUE, reversed = TRUE, details = TRUE)
#save the results to a file. 
write.csv(as.data.frame(pairwiseletters$emmeans), file = "OptimizationALLPairwiseEmmeans.csv")
write.csv(as.data.frame(pairwiseletters$comparisons), file = "OptimizationALLPairwiseComparisons.csv")
#save the workspace
save.image("OptimizationAfterPairwise.RData")

#MSP1-19, 100 ug/mL, 
colnames(AHHHH.df[17])
fullmodel <- lmer(X19_1.PfMSP1.19.100ug.ml_1 ~ slide_type + blocking_buffer 
                  + block_dilution + print_buffer + (1|sample), 
                  REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#plot residuals - they look mostly good, some heteroskedasticity though and data looks linear
plot(fitted(fullmodel),residuals(fullmodel))
#check normality - looks good
hist(residuals(fullmodel))

#prepare a better organized summary of the model and export to a file.
#summary can only be written to a text file, and doesn't keep columns organized
write.csv(tidy(fullmodel), file = "MSP1-19.100.LMERTidy.csv")
write.csv(augment(fullmodel), file = "MSP1-19.100.LMERAug.csv")

#test for main effects and interactions with Likelihood Ratio Test?
intmodel <- lmer(X19_1.PfMSP1.19.100ug.ml_1 ~ slide_type * blocking_buffer 
                 * block_dilution * print_buffer + (1|sample), 
                 REML = FALSE, data = AHHHH.df)
summary(intmodel)

#Plot residuals 
png(filename = "MSP1.19.100.Residuals.tif", width = 7.8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

plot(fitted(intmodel),residuals(intmodel),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
abline(a=0, b=0)
hist(residuals(intmodel), xlab = "Residuals")

graphics.off()

r.squaredGLMM(intmodel)

write.csv(tidy(intmodel), file = "MSP1.19.100.LMERTidyINT.csv")
write.csv(augment(intmodel), file = "MSP1.19.100.LMERAugINT.csv")

anova(fullmodel, intmodel)

anova(intmodel)

#get significance and do pairwise comparisons using emmeans
#default adjustment for multiple comparisons is Tukey
#main effects first, however we are not interested in main effects :(
emmeans(fullmodel, pairwise ~ print_buffer)
emmeans(fullmodel, pairwise ~ slide_type)
emmeans(fullmodel, pairwise ~ blocking_buffer)
emmeans(fullmodel, pairwise ~ block_dilution)

#plot the data
png(filename = "MSP1-19.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type)

graphics.off()

#get letters for all pairwise comparisons, interaction model
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#this took several hours (>4) to run, and eventually finished. there are too many 
#combinations (480 pairwise), but it worked!!! 
#the contrast function used by cld function automatically uses Tukey adjustment
#for multiple comparisons
pairwiseletters <- cld(emm.intmodel, by = NULL, Letters = LETTERS, sort = TRUE, reversed = TRUE, details = TRUE)
#save the results to a file. 
write.csv(as.data.frame(pairwiseletters$emmeans), file = "MSP1-19.100.ALLPairwiseEmmeans.csv")
write.csv(as.data.frame(pairwiseletters$comparisons), file = "MSP1-19.100.ALLPairwiseComparisons.csv")
#above, forgot to change the name of pairwiseletters to be antigen specific!
#so this pairwise overwrote the previous one for AMA1

#save the workspace
save.image("OptimizationAfterPairwiseMSP1-19.RData")

#Hyp2, 100 ug/mL, 
colnames(AHHHH.df[23])
fullmodel <- lmer(X25_1.Hyp2.100ug.ml_1 ~ slide_type + blocking_buffer 
                  + block_dilution + print_buffer + (1|sample), 
                  REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#plot residuals - they look mostly good, some heteroskedasticity though and data looks linear
plot(fitted(fullmodel),residuals(fullmodel))
#check normality - looks good
hist(residuals(fullmodel))

#prepare a better organized summary of the model and export to a file.
#summary can only be written to a text file, and doesn't keep columns organized
write.csv(tidy(fullmodel), file = "Hyp2.100.LMERTidy.csv")
write.csv(augment(fullmodel), file = "Hyp2.100.LMERAug.csv")

#test for main effects and interactions with Likelihood Ratio Test?
intmodel <- lmer(X25_1.Hyp2.100ug.ml_1 ~ slide_type * blocking_buffer 
                 * block_dilution * print_buffer + (1|sample), 
                 REML = FALSE, data = AHHHH.df)
summary(intmodel)

#Plot residuals 
png(filename = "Hyp2.Residuals.tif", width = 7.8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

plot(fitted(intmodel),residuals(intmodel),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
abline(a=0, b=0)
hist(residuals(intmodel), xlab = "Residuals")

graphics.off()

r.squaredGLMM(intmodel)

write.csv(tidy(intmodel), file = "Hyp2.100.LMERTidyINT.csv")
write.csv(augment(intmodel), file = "Hyp2.100.LMERAugINT.csv")

anova(fullmodel, intmodel)

anova(intmodel)

#plot the model data
png(filename = "Hyp2.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type)

graphics.off()

#get letters for all pairwise comparisons, interaction model
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#this took several hours (>6.5) to run, and eventually finished. there are too many 
#combinations (480 pairwise), but it worked!!! 
#the contrast function used by cld function automatically uses Tukey adjustment
#for multiple comparisons
Hyp2pairwiseletters <- cld(emm.intmodel, by = NULL, Letters = LETTERS, sort = TRUE, reversed = TRUE, details = TRUE)
#save the results to a file. 
write.csv(as.data.frame(Hyp2pairwiseletters$emmeans), file = "Hyp2.100.ALLPairwiseEmmeans.csv")
write.csv(as.data.frame(Hyp2pairwiseletters$comparisons), file = "Hyp2.100.ALLPairwiseComparisons.csv")


######### Plotting the data from the successful combos in pairwiseletters #######

#import the data which I decided was good in excel
AMA1.100.best <- read.csv("AMA1.100.BEST.csv")
MSP1.19.100.best <- read.csv("MSP1-19.100.BEST.csv")

#add rowid as a column to get a number for each combination
AMA1.100.best <- tibble::rowid_to_column(AMA1.100.best)
MSP1.19.100.best <- tibble::rowid_to_column(MSP1.19.100.best)

#merge with the AHHHH.df to get original values for all three samples
#some are NA! I didn't realize that samples with NA would be included
AMA1.100.exp <- merge(AMA1.100.best, AHHHH.df, all.x = TRUE, sort = FALSE)
MSP1.19.100.exp <- merge(MSP1.19.100.best, AHHHH.df, all.x = TRUE, sort = FALSE)

#keep same order that we have in the excel file (sorted highest to lowest EMM)
#maybe try to use rownumber as factor, sorted by EMM? and then put in better labels later?
rowidorder <- c(as.character(AMA1.100.exp$rowid))
AMA1.100.exp$rowid <- factor(AMA1.100.exp$rowid, levels = rev(unique(rowidorder)))

rowidorderMSP <- c(as.character(MSP1.19.100.exp$rowid))
MSP1.19.100.exp$rowid <- factor(MSP1.19.100.exp$rowid, levels = rev(unique(rowidorderMSP)))

#COORD_FLIP!!! Plot every point for each combo, different colors for each sample, 
#AMA1 100
png(filename = paste0("AMA1.100.best.dot.tif"), width = 3, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AMA1.100.exp, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + coord_flip() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, hjust=1, size = 10))

graphics.off()

#MSP1.19 100
png(filename = paste0("MSP1_19.100.best.dot.tif"), width = 3, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(MSP1.19.100.exp, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + coord_flip() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, hjust=1, size = 10))

graphics.off()

#Same plots - NOT COORD_FLIP to put on PPT slides 
#AMA1 100
AMA1.100.exp$rowid <- factor(AMA1.100.exp$rowid, levels = unique(rowidorder))

png(filename = paste0("H.AMA1.100.best.dot.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AMA1.100.exp, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

#MSP1.19 100
MSP1.19.100.exp$rowid <- factor(MSP1.19.100.exp$rowid, levels = unique(rowidorderMSP))

png(filename = paste0("H.MSP1_19.100.best.dot.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(MSP1.19.100.exp, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text(vjust = 1, size = 10))

graphics.off()

#Which combinations are present in multiple antigens? 
AMA1.100.Combo <- AMA1.100.best[,2:5]
MSP1.19.100.Combo <- MSP1.19.100.best[,2:5]

ComboMult <- merge(AMA1.100.Combo, MSP1.19.100.Combo, sort = FALSE)
#12 combos are present in AMA1.100 and MSP1.19.100 best lists